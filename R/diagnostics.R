################################################################################
## Diagnostics
##
## panelcoder() records every estimation problem it can detect as one row of a
## `diagnostics` data frame (columns type, severity, message, parameter, value)
## and summarizes them in a one-row `status` data frame. Both are added to the
## pcOutput as named elements, so the positional elements 1-8 (including the
## legacy single-string warning [[7]] and error [[8]]) are unchanged.
##
## Severities:
##   error          estimation failed or did not converge; results unusable
##   inadmissible   a non-positive definite latent or residual covariance matrix
##   local_optimum  the best solution was not replicated across random starts
##   se             standard errors are missing, zero, or flagged as untrustworthy
##   boundary       a constraint panelCodeR imposes (limits = TRUE) is active
##   info           anything else worth knowing
################################################################################

## Tolerances for active-bound detection: a variance is "at 0" when it is
## smaller than .boundVarTol times that variable's total observed variance; a
## correlation is "at the .99 bound" when |r| >= .99 - .boundCorTol.
.boundVarTol <- 1e-4
.boundCorTol <- 1e-3

## Random-starts replication: solutions count as replicating the best one when
## within this distance of it (Mplus fit function values are printed to four
## decimals; loglikelihood values to three).
.startsFitTol <- 1e-3
.startsLLTol <- .01

## Status flags that each diagnostic type sets to FALSE
.notConvergedTypes <- c("no_convergence", "not_terminated", "mplus_input_error",
                        "sparse", "r_error", "no_estimates")
.seUntrustedTypes <- c("se_untrusted", "se_failed", "se_undefined")
.inadmissibleTypes <- c("psi_npd", "theta_npd", "post_check")

.emptyDiagnostics <- function() {
    data.frame(type = character(0),
               severity = character(0),
               message = character(0),
               parameter = character(0),
               value = numeric(0),
               stringsAsFactors = FALSE)
}

.diagRow <- function(type, severity, message, parameter = NA_character_,
                     value = NA_real_) {
    data.frame(type = type,
               severity = severity,
               message = message,
               parameter = as.character(parameter),
               value = as.numeric(value),
               stringsAsFactors = FALSE)
}

.bindDiagnostics <- function(...) {
    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, list(...))
    if (length(parts) == 0) return(.emptyDiagnostics())
    out <- do.call(rbind, parts)
    out <- out[!duplicated(out[, c("type", "message", "parameter")]), ,
               drop = FALSE]
    rownames(out) <- NULL
    out
}

.normalizeWs <- function(x) trimws(gsub("\\s+", " ", x))

## Convert Mplus Fortran-style numbers ("0.647D-12") to numeric
.mplusNum <- function(x) {
    suppressWarnings(as.numeric(sub("D", "E", x, fixed = TRUE)))
}

## Split Mplus output into paragraphs (runs of non-blank lines). Returns the
## whitespace-normalized text of each paragraph.
.mplusParagraphs <- function(lines) {
    if (length(lines) == 0) return(character(0))
    blank <- !grepl("\\S", lines)
    id <- cumsum(blank)
    keep <- !blank
    if (!any(keep)) return(character(0))
    paras <- tapply(lines[keep], id[keep], paste, collapse = " ")
    unname(.normalizeWs(as.character(paras)))
}

## Mplus condition number in a paragraph ("THE CONDITION NUMBER IS -0.674D-05.")
.extractCondNumber <- function(text) {
    m <- regmatches(text, regexpr("CONDITION NUMBER IS\\s+-?[0-9.]+D[-+][0-9]+",
                                  text))
    if (length(m) == 0) return(NA_real_)
    .mplusNum(sub("CONDITION NUMBER IS\\s+", "", m))
}

## "PROBLEM INVOLVING THE FOLLOWING PARAMETER: Parameter 6, X_2"
.extractProblemParameter <- function(text) {
    m <- regmatches(text, regexpr(
        "PROBLEM INVOLVING THE FOLLOWING PARAMETER:\\s*Parameter [0-9]+, [^ ]+( [A-Z_0-9]+ [A-Z_0-9]+)?",
        text))
    if (length(m) == 0) {
        m <- regmatches(text, regexpr(
            "PROBLEM INVOLVING PARAMETER\\s+[0-9]+", text))
        if (length(m) == 0) return(NA_character_)
        return(paste("Parameter", sub("PROBLEM INVOLVING PARAMETER\\s+", "", m)))
    }
    m <- sub("PROBLEM INVOLVING THE FOLLOWING PARAMETER:\\s*", "", m)
    ## Drop any following sentence that was pulled in
    sub("\\s+(THE|THIS)\\s.*$", "", sub("\\.$", "", m))
}

## "PROBLEM INVOLVING VARIABLE S_X_1."
.extractProblemVariable <- function(text) {
    m <- regmatches(text, regexpr("PROBLEM INVOLVING VARIABLE\\s+[^ .]+", text))
    if (length(m) == 0) return(NA_character_)
    sub("PROBLEM INVOLVING VARIABLE\\s+", "", m)
}

## Classify a single (normalized) Mplus paragraph. Returns a diagnostics data
## frame with one row per known message in the paragraph (zero rows when none).
.classifyMplusParagraph <- function(p, nextP = "") {
    out <- list()
    add <- function(row) out[[length(out) + 1]] <<- row
    if (grepl("^\\*\\*\\* (FATAL )?ERROR", p)) {
        add(.diagRow("mplus_input_error", "error", p))
    }
    if (grepl("MAY NOT BE TRUSTWORTHY FOR SOME PARAMETERS DUE TO A NON-POSITIVE DEFINITE FIRST-ORDER DERIVATIVE PRODUCT MATRIX",
              p, fixed = TRUE)) {
        cn <- .extractCondNumber(p)
        if (is.na(cn)) cn <- .extractCondNumber(nextP)
        add(.diagRow("se_untrusted", "se",
                        "Standard errors may not be trustworthy (non-positive definite first-order derivative product matrix)",
                        parameter = .extractProblemParameter(paste(p, nextP)),
                        value = cn))
    }
    if (grepl("STANDARD ERRORS COULD NOT BE COMPUTED\\. THE MLF STANDARD ERRORS WERE COMPUTED INSTEAD", p)) {
        cn <- .extractCondNumber(p)
        add(.diagRow("se_untrusted", "se",
                        "ML/MLR standard errors could not be computed; MLF standard errors were used instead",
                        parameter = .extractProblemParameter(p),
                        value = cn))
    }
    if (grepl("THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES COULD NOT BE COMPUTED",
              p, fixed = TRUE)) {
        cn <- .extractCondNumber(p)
        if (is.na(cn)) cn <- .extractCondNumber(nextP)
        add(.diagRow("se_failed", "error",
                        "Standard errors could not be computed (the model may not be identified)",
                        parameter = .extractProblemParameter(paste(p, nextP)),
                        value = cn))
    }
    if (grepl("THE LATENT VARIABLE COVARIANCE MATRIX (PSI) IS NOT POSITIVE DEFINITE",
              p, fixed = TRUE)) {
        v <- .extractProblemVariable(p)
        msg <- "Latent variable covariance matrix (PSI) is not positive definite"
        if (!is.na(v)) msg <- paste0(msg, " (problem involving ", v, ")")
        add(.diagRow("psi_npd", "inadmissible", msg, parameter = v))
    }
    if (grepl("THE RESIDUAL COVARIANCE MATRIX (THETA) IS NOT POSITIVE DEFINITE",
              p, fixed = TRUE)) {
        v <- .extractProblemVariable(p)
        msg <- "Residual covariance matrix (THETA) is not positive definite"
        if (!is.na(v)) msg <- paste0(msg, " (problem involving ", v, ")")
        add(.diagRow("theta_npd", "inadmissible", msg, parameter = v))
    }
    if (grepl("THE BEST LOGLIKELIHOOD VALUE WAS NOT REPLICATED", p, fixed = TRUE)) {
        add(.diagRow("ll_not_replicated", "local_optimum",
                        "The best loglikelihood value was not replicated; the solution may be a local maximum"))
    }
    if (grepl("THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY", p, fixed = TRUE)) {
        reason <- sub("^.*THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY\\.?\\s*", "", p)
        reason <- sub("^DUE TO\\s+", "", reason)
        reason <- sub("\\..*$", "", reason)
        msg <- "Model estimation did not terminate normally"
        if (nzchar(reason)) msg <- paste0(msg, " due to ", tolower(reason))
        add(.diagRow("not_terminated", "error", msg,
                        value = .extractCondNumber(p)))
    }
    if (grepl("NO CONVERGENCE\\.", p)) {
        variant <- sub("^.*?NO CONVERGENCE\\.\\s*", "", p, perl = TRUE)
        variant <- sub("\\..*$", "", variant)
        msg <- "No convergence"
        if (nzchar(variant)) msg <- paste0(msg, ": ", tolower(variant))
        add(.diagRow("no_convergence", "error", msg))
    }
    if (grepl("THE MISSING DATA EM ALGORITHM FOR THE H1 MODEL HAS NOT CONVERGED",
              p, fixed = TRUE)) {
        add(.diagRow("sparse", "error",
                        "The missing data EM algorithm for the H1 model has not converged (sparse data)"))
    }
    if (grepl("THE COVARIANCE COVERAGE FALLS BELOW THE SPECIFIED LIMIT",
              p, fixed = TRUE)) {
        add(.diagRow("coverage", "info",
                        "Covariance coverage falls below the specified limit"))
    }
    do.call(.bindDiagnostics, out)
}

## Random-starts table. Returns list(found, n, nReplicated, best, values)
## or NULL when the output has no table. Fit function values are minimized;
## loglikelihood values are maximized.
.mplusStartsTable <- function(lines) {
    hdr <- grep("(Fit function|[Ll]oglikelihood) values at local maxima", lines)
    if (length(hdr) == 0) return(NULL)
    hdr <- hdr[length(hdr)]
    isLL <- grepl("oglikelihood", lines[hdr])
    vals <- numeric(0)
    i <- hdr + 1
    started <- FALSE
    while (i <= length(lines)) {
        ln <- lines[i]
        if (!grepl("\\S", ln)) {
            if (started) break
            i <- i + 1
            next
        }
        m <- regmatches(ln, regexpr("^\\s*-?[0-9]+\\.[0-9]+(D[-+][0-9]+)?(?=\\s)",
                                    ln, perl = TRUE))
        if (length(m) == 0) break
        vals <- c(vals, .mplusNum(trimws(m)))
        started <- TRUE
        i <- i + 1
    }
    vals <- vals[is.finite(vals)]
    ## Starts that did not converge ("Unperturbed starting value run did not
    ## converge." and "4 perturbed starting value run(s) did not converge.")
    after <- lines[seq.int(hdr, length(lines))]
    nFailed <- sum(grepl("^\\s*Unperturbed starting value run did not converge",
                         after))
    m <- regmatches(after, regexpr(
        "^\\s*[0-9]+ (perturbed )?starting value run\\(s\\) did not converge",
        after))
    if (length(m) > 0) nFailed <- nFailed + sum(as.numeric(sub("^\\s*([0-9]+).*", "\\1", m)))
    if (length(vals) == 0) {
        return(list(n = 0, nReplicated = 0, best = NA_real_, isLL = isLL,
                    nFailed = nFailed))
    }
    best <- if (isLL) max(vals) else min(vals)
    tol <- if (isLL) .startsLLTol else .startsFitTol
    list(n = length(vals),
         nReplicated = sum(abs(vals - best) <= tol + 1e-12),
         best = best,
         isLL = isLL,
         nFailed = nFailed)
}

## Was the best solution replicated across random starts? NA when random
## starts were not used or fewer than two starts converged (replication
## cannot be checked).
.mplusLLReplicated <- function(lines) {
    text <- paste(.mplusParagraphs(lines), collapse = " ")
    if (grepl("THE BEST LOGLIKELIHOOD VALUE WAS NOT REPLICATED", text, fixed = TRUE)) {
        return(FALSE)
    }
    if (grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", text, fixed = TRUE)) {
        return(TRUE)
    }
    tab <- .mplusStartsTable(lines)
    if (is.null(tab) || tab$n < 2) return(NA)
    tab$nReplicated >= 2
}

## Smallest covariance coverage value in the "PROPORTION OF DATA PRESENT"
## table (NA if not found)
.mplusMinCoverage <- function(lines) {
    start <- grep("PROPORTION OF DATA PRESENT", lines, fixed = TRUE)
    if (length(start) == 0) return(NA_real_)
    i <- start[1] + 1
    vals <- numeric(0)
    while (i <= length(lines) && !grepl("^[A-Z]", lines[i])) {
        ln <- lines[i]
        if (grepl("^\\s+[A-Z][A-Z0-9_]*\\s+[0-9]", ln)) {
            nums <- regmatches(ln, gregexpr("(?<=\\s)[0-9]+\\.[0-9]+", ln,
                                            perl = TRUE))[[1]]
            vals <- c(vals, as.numeric(nums))
        }
        i <- i + 1
    }
    if (length(vals) == 0) NA_real_ else min(vals)
}

## Parameter estimates and SEs from an Mplus fit, matched to the labels in
## the panelCodeR parameter table. Only free, labeled ~, ~~, and =~ rows.
.mplusLabelEstimates <- function(model, fit) {
    est <- fit$results$parameters$unstandardized
    if (is.null(est) || is.null(model)) return(NULL)
    rows <- which(model$free != 0 & model$label != "" &
                  model$op %in% c("~", "~~", "=~"))
    vals <- rep(NA_real_, length(rows))
    ses <- rep(NA_real_, length(rows))
    hdr <- toupper(est$paramHeader)
    prm <- toupper(est$param)
    for (k in seq_along(rows)) {
        r <- rows[k]
        lhs <- toupper(gsub(".", "_", model$lhs[r], fixed = TRUE))
        rhs <- toupper(gsub(".", "_", model$rhs[r], fixed = TRUE))
        op <- model$op[r]
        if (op == "~") {
            hit <- which(hdr == paste0(lhs, ".ON") & prm == rhs)
        } else if (op == "=~") {
            hit <- which(hdr == paste0(lhs, ".BY") & prm == rhs)
        } else if (lhs == rhs) {
            hit <- which(hdr %in% c("VARIANCES", "RESIDUAL.VARIANCES") &
                         prm == lhs)
        } else {
            hit <- which((hdr == paste0(lhs, ".WITH") & prm == rhs) |
                         (hdr == paste0(rhs, ".WITH") & prm == lhs))
        }
        if (length(hit) > 0) {
            vals[k] <- est$est[hit[1]]
            if (!is.null(est$se)) ses[k] <- est$se[hit[1]]
        }
    }
    data.frame(label = model$label[rows],
               lhs = model$lhs[rows],
               op = model$op[rows],
               rhs = model$rhs[rows],
               est = vals,
               se = ses,
               stringsAsFactors = FALSE)
}

## Parameter estimates and SEs from a lavaan fit, by label
.lavaanLabelEstimates <- function(model, fit) {
    pt <- tryCatch(lavaan::parTable(fit), error = function(e) NULL)
    if (is.null(pt)) return(NULL)
    keep <- pt$label != "" & pt$free > 0 & pt$op %in% c("~", "~~", "=~")
    data.frame(label = pt$label[keep],
               lhs = pt$lhs[keep],
               op = pt$op[keep],
               rhs = pt$rhs[keep],
               est = pt$est[keep],
               se = if (is.null(pt$se)) NA_real_ else pt$se[keep],
               stringsAsFactors = FALSE)
}

## Free parameters whose SE is missing, 0, or 999 (Mplus's marker for an
## undefined value), excluding parameters that depend on others through `==`.
.undefinedSEs <- function(labelEst, model) {
    if (is.null(labelEst) || nrow(labelEst) == 0) return(.emptyDiagnostics())
    deps <- if (is.null(model)) character(0) else .dependentLabels(model)
    le <- labelEst[!(labelEst$label %in% deps) & !is.na(labelEst$est), ,
                   drop = FALSE]
    bad <- is.na(le$se) | le$se == 0 | le$se == 999
    if (!any(bad)) return(.emptyDiagnostics())
    le <- le[bad, , drop = FALSE]
    msgs <- paste0("Standard error of ", le$label, " (", le$lhs, " ", le$op,
                   " ", le$rhs, ") is ",
                   ifelse(is.na(le$se), "missing",
                          ifelse(le$se == 999, "undefined (999)", "0")))
    .bindDiagnostics(.diagRow("se_undefined", "se", msgs, parameter = le$label,
                              value = le$se))
}

## Readable name of the parameter with a given label
.describeLabel <- function(label, model, info, roles) {
    r <- which(model$label == label & model$op %in% c("~", "~~", "=~"))
    if (length(r) == 0) return(label)
    r <- r[1]
    role <- roles[r]
    if (is.na(role)) return(label)
    var <- sub("^.*_(x|y|xy|yx)$", "\\1", role)
    base <- sub("_(x|y|xy|yx)$", "", role)
    words <- switch(base,
        trait = "stable-trait variance",
        ar = "autoregressive variance",
        innov = "innovation variance",
        state = "state variance",
        slope = "slope variance",
        resid = "residual variance",
        cov_trait = "stable-trait correlation",
        cov_ar = "autoregressive correlation",
        cov_innov = "innovation correlation",
        cov_state = "state correlation",
        cov_slope = "slope correlation",
        NULL)
    if (is.null(words)) return(label)
    if (var %in% c("x", "y") && !is.null(info[[var]]$name)) {
        words <- paste0(words, " of ", info[[var]]$name)
        ## Add the wave when this parameter is not shared across waves
        wave <- .parsePcName(model$lhs[r], info)$wave
        shared <- any(model$op == "==" & model$rhs == label)
        if (!is.na(wave) && !shared) {
            words <- paste0(words, " (wave ", wave, ")")
        }
    }
    words
}

## Total observed variance of the variable a label belongs to
.labelScale <- function(label, model, info, data) {
    if (is.null(data)) return(1)
    allVar <- mean(vapply(data, function(v) {
        if (is.numeric(v)) stats::var(v, na.rm = TRUE) else NA_real_
    }, numeric(1)), na.rm = TRUE)
    if (!is.finite(allVar) || allVar <= 0) allVar <- 1
    r <- which(model$label == label & model$op == "~~")
    if (length(r) == 0) return(allVar)
    var <- .parsePcName(model$lhs[r[1]], info)$var
    if (is.na(var) || is.null(info[[var]]$name)) return(allVar)
    stem <- info[[var]]$name
    cols <- names(data)[startsWith(names(data), paste0(stem, "_"))]
    if (length(cols) == 0) return(allVar)
    s <- mean(vapply(data[cols], stats::var, numeric(1), na.rm = TRUE),
              na.rm = TRUE)
    if (!is.finite(s) || s <= 0) allVar else s
}

## Inequality constraints (from limits = TRUE) that the solution sits on
.activeBounds <- function(model, labelEst, info, data = NULL) {
    if (is.null(model) || is.null(labelEst) || nrow(labelEst) == 0) {
        return(.emptyDiagnostics())
    }
    cons <- which(model$op %in% c(">", "<"))
    if (length(cons) == 0) return(.emptyDiagnostics())
    ests <- labelEst$est[!is.na(labelEst$est)]
    names(ests) <- labelEst$label[!is.na(labelEst$est)]
    ests <- ests[!duplicated(names(ests))]
    env <- list2env(as.list(ests), parent = baseenv())
    evalSide <- function(expr) {
        val <- tryCatch(eval(parse(text = expr), envir = env),
                        error = function(e) NA_real_)
        if (!is.numeric(val) || length(val) != 1) NA_real_ else val
    }
    roles <- .startValRoles(model, info)
    out <- list()
    for (r in cons) {
        lhsExpr <- model$lhs[r]
        rhsExpr <- as.character(model$rhs[r])
        lhsVal <- evalSide(lhsExpr)
        rhsVal <- evalSide(rhsExpr)
        if (is.na(lhsVal) || is.na(rhsVal)) next
        sqrtTerms <- regmatches(rhsExpr,
                                gregexpr("sqrt\\(([A-Za-z0-9_.]+)\\)", rhsExpr))[[1]]
        if (length(sqrtTerms) == 2) {
            ## Correlation bound: lhs (a covariance) relative to sd * sd
            vars <- sub("sqrt\\((.*)\\)", "\\1", sqrtTerms)
            sds <- vapply(vars, function(v) evalSide(paste0("sqrt(", v, ")")),
                          numeric(1))
            if (any(is.na(sds)) || any(sds <= 0)) next
            corr <- lhsVal / prod(sds)
            k <- abs(rhsVal) / prod(sds)
            if (abs(corr) >= k - .boundCorTol &&
                sign(corr) == sign(rhsVal)) {
                what <- .describeLabel(lhsExpr, model, info, roles)
                if (identical(what, lhsExpr)) what <- paste("correlation", lhsExpr)
                bound <- sub("^0", "", format(round(k, 2), nsmall = 2))
                if (rhsVal < 0) bound <- paste0("-", bound)
                out[[length(out) + 1]] <- .diagRow(
                    "boundary", "boundary",
                    paste0(what, " at the ", bound, " bound"),
                    parameter = lhsExpr, value = corr)
            }
        } else if (rhsVal == 0) {
            scale <- .labelScale(lhsExpr, model, info, data)
            if (abs(lhsVal) < .boundVarTol * scale) {
                what <- .describeLabel(lhsExpr, model, info, roles)
                out[[length(out) + 1]] <- .diagRow(
                    "boundary", "boundary", paste0(what, " at 0"),
                    parameter = lhsExpr, value = lhsVal)
            }
        } else {
            if (abs(lhsVal - rhsVal) <= .boundCorTol * max(abs(rhsVal), 1e-8)) {
                out[[length(out) + 1]] <- .diagRow(
                    "boundary", "boundary",
                    paste0("constraint ", lhsExpr, " ", model$op[r], " ",
                           rhsExpr, " is active"),
                    parameter = lhsExpr, value = lhsVal)
            }
        }
    }
    do.call(.bindDiagnostics, out)
}

## All diagnostics that can be read from Mplus output text (and, when given,
## the MplusAutomation fit object and the panelCodeR parameter table).
.mplusDiagnostics <- function(outputLines, fit = NULL, model = NULL) {
    out <- list()
    paras <- .mplusParagraphs(outputLines)
    for (i in seq_along(paras)) {
        nextP <- if (i < length(paras)) paras[i + 1] else ""
        out[[length(out) + 1]] <- .classifyMplusParagraph(paras[i], nextP)
    }

    ## Random starts: local maxima
    out <- list(do.call(.bindDiagnostics, out))
    if (!("ll_not_replicated" %in% out[[1]]$type)) {
        tab <- .mplusStartsTable(outputLines)
        if (!is.null(tab) && tab$n == 1) {
            ## Only one start converged: nothing to replicate against
            out[[length(out) + 1]] <- .diagRow(
                "starts_insufficient", "info",
                paste0("Only 1 of ", 1 + tab$nFailed,
                       " starts converged; replication of the best solution could not be checked"),
                value = 1)
        } else if (!is.null(tab) && tab$n >= 2 && tab$nReplicated < 2) {
            what <- if (tab$isLL) "loglikelihood" else "fit function value"
            out[[length(out) + 1]] <- .diagRow(
                "ll_not_replicated", "local_optimum",
                paste0("The best ", what, " was not replicated across random starts (",
                       tab$n, " converged starts)"),
                value = tab$nReplicated)
        }
    }

    ## Coverage value
    hasCov <- out[[1]]$type == "coverage"
    if (any(hasCov)) {
        out[[1]]$value[hasCov] <- .mplusMinCoverage(outputLines)
    }

    ## Messages MplusAutomation collected that are not covered above
    if (!is.null(fit$results)) {
        msgs <- c(fit$results$warnings, fit$results$errors)
        for (m in msgs) {
            p <- .normalizeWs(paste(unlist(m), collapse = " "))
            if (!nzchar(p)) next
            known <- out[[1]]$message[out[[1]]$type == "mplus_input_error"]
            if (nrow(.classifyMplusParagraph(p)) == 0 &&
                !any(vapply(known, grepl, logical(1), pattern = p,
                            fixed = TRUE))) {
                out[[length(out) + 1]] <- .diagRow("mplus_warning", "info", p)
            }
        }
    }

    ## Free parameters with SE of 0 or 999 (only for converged solutions with
    ## SEs, otherwise every parameter would be listed)
    skipSE <- any(out[[1]]$type %in% c("se_failed", .notConvergedTypes))
    if (!skipSE && !is.null(model)) {
        out[[length(out) + 1]] <- .undefinedSEs(.mplusLabelEstimates(model, fit),
                                                model)
    }
    do.call(.bindDiagnostics, out)
}

## lavaan: recorded R warnings plus convergence, admissibility and SE checks
.lavaanDiagnostics <- function(fit, rWarnings = character(0), model = NULL) {
    out <- list()
    for (w in unique(rWarnings)) {
        out[[length(out) + 1]] <- .diagRow("r_warning", "info", .normalizeWs(w))
    }
    if (is.null(fit)) return(do.call(.bindDiagnostics, out))

    conv <- tryCatch(lavaan::lavInspect(fit, "converged"),
                     error = function(e) NA)
    if (isFALSE(conv)) {
        out[[length(out) + 1]] <- .diagRow("no_convergence", "error",
                                           "lavaan did not converge")
    }
    if (isTRUE(conv)) {
        post <- tryCatch(suppressWarnings(lavaan::lavInspect(fit, "post.check")),
                         error = function(e) NA)
        if (isFALSE(post)) {
            why <- rWarnings[grepl("negative|not positive definite", rWarnings) &
                             !grepl("vcov", rWarnings)]
            msg <- "Inadmissible solution (lavaan post.check failed: negative variance or non-positive definite matrix)"
            if (length(why) > 0) msg <- paste0(msg, ": ", .normalizeWs(why[1]))
            out[[length(out) + 1]] <- .diagRow("post_check", "inadmissible", msg)
        }
    }
    vcovWarn <- rWarnings[grepl("vcov", rWarnings) &
                          grepl("positive definite", rWarnings)]
    if (length(vcovWarn) > 0) {
        ev <- regmatches(vcovWarn[1],
                         regexpr("-?[0-9.]+e[-+][0-9]+", vcovWarn[1]))
        out[[length(out) + 1]] <- .diagRow(
            "se_untrusted", "se",
            "Variance-covariance matrix of the estimates is not positive definite",
            value = if (length(ev) > 0) as.numeric(ev) else NA_real_)
    }
    seWarn <- any(grepl("could not compute standard errors", rWarnings))
    le <- .lavaanLabelEstimates(model, fit)
    deps <- if (is.null(model)) character(0) else .dependentLabels(model)
    indep <- if (is.null(le)) NULL else le[!(le$label %in% deps), , drop = FALSE]
    allMissing <- !is.null(indep) && nrow(indep) > 0 && all(!is.finite(indep$se))
    if (isTRUE(conv) && (seWarn || allMissing)) {
        out[[length(out) + 1]] <- .diagRow("se_failed", "error",
                                           "Standard errors could not be computed")
    } else if (isTRUE(conv)) {
        out[[length(out) + 1]] <- .undefinedSEs(le, model)
    }
    do.call(.bindDiagnostics, out)
}

## One-row summary of the diagnostics
.diagStatus <- function(diag, hasEstimates, llReplicated = NA,
                        minCoverage = NA_real_) {
    types <- diag$type
    condRows <- diag$type %in% c("se_untrusted", "se_failed", "not_terminated") &
        !is.na(diag$value)
    data.frame(
        converged = hasEstimates && !any(types %in% .notConvergedTypes),
        seTrusted = if (hasEstimates) !any(types %in% .seUntrustedTypes) else NA,
        admissible = if (hasEstimates) !any(types %in% .inadmissibleTypes) else NA,
        llReplicated = as.logical(llReplicated),
        boundary = any(types == "boundary"),
        condNumber = if (any(condRows)) min(diag$value[condRows]) else NA_real_,
        minCoverage = as.numeric(minCoverage),
        stringsAsFactors = FALSE)
}

## Smallest lavaan coverage value (NA when not available)
.lavaanMinCoverage <- function(fit) {
    cv <- tryCatch(suppressWarnings(lavaan::lavInspect(fit, "coverage")),
                   error = function(e) NULL)
    if (is.null(cv) || length(cv) == 0) return(NA_real_)
    if (is.list(cv)) cv <- unlist(cv)
    suppressWarnings(min(cv, na.rm = TRUE))
}

## Legacy single warning / error strings for pcOutput[[7]] and [[8]]. These
## reproduce the messages panelcoder() gave before diagnostics were added, so
## code that reads [[7]] and [[8]] sees the same values.
.legacyMplusMessages <- function(diag, lines, success = TRUE) {
    types <- diag$type
    warningM <- NULL
    errorM <- NULL
    inputError <- .mplusErrorFromLines(lines)
    if (!is.null(inputError)) {
        return(list(warningM = NULL, errorM = inputError, success = FALSE))
    }
    if (any(types == "se_untrusted" &
            grepl("first-order derivative", diag$message))) {
        warningM <- "Non-positive definite matrix"
    }
    if ("se_failed" %in% types) {
        success <- FALSE
        errorM <- "Standard errors could not be computed"
    }
    if ("no_convergence" %in% types) {
        success <- FALSE
        errorM <- "No convergence"
    }
    if ("sparse" %in% types) {
        success <- FALSE
        errorM <- "No convergence; Sparse Matrix"
    }
    if ("coverage" %in% types) {
        warningM <- "Coverage warning"
    }
    list(warningM = warningM, errorM = errorM, success = success)
}

## Run the Mplus checks on a finished (or failed) Mplus run and return the
## pieces panelcoder() needs: success flag, legacy warning/error strings,
## diagnostics and status.
.mplusRunResult <- function(fit, file_out, rWarnings = character(0),
                            runError = NULL, model = NULL, info = NULL,
                            data = NULL) {
    lines <- fit$results$output
    if (is.null(lines) && !is.null(file_out) && file.exists(file_out)) {
        lines <- readLines(file_out, warn = FALSE)
    }
    rRows <- lapply(unique(rWarnings), function(w) {
        .diagRow("r_warning", "info", .normalizeWs(w))
    })
    warningM <- if (length(rWarnings) > 0) rWarnings[length(rWarnings)] else NULL

    if (is.null(lines)) {
        msg <- if (!is.null(runError)) runError else check_mplus_error(file_out)
        diag <- do.call(.bindDiagnostics, c(
            rRows,
            list(.diagRow(if (!is.null(runError)) "r_error" else "mplus_input_error",
                          "error", msg))))
        return(list(fit = fit, success = FALSE, warningM = warningM,
                    errorM = msg, diagnostics = diag,
                    status = .diagStatus(diag, FALSE)))
    }

    mDiag <- .mplusDiagnostics(lines, fit, model)
    legacy <- .legacyMplusMessages(mDiag, lines)
    success <- legacy$success
    errorM <- legacy$errorM
    if (!is.null(legacy$warningM)) warningM <- legacy$warningM

    extra <- list()
    if (!is.null(runError)) {
        success <- FALSE
        errorM <- runError
        extra[[1]] <- .diagRow("r_error", "error", .normalizeWs(runError))
    }

    hasEstimates <- !is.null(fit$results$parameters$unstandardized)
    if (success && !hasEstimates) {
        success <- FALSE
        mplusError <- check_mplus_error(file_out)
        errorM <- if (!is.null(mplusError)) {
            mplusError
        } else {
            paste0("Mplus did not produce parameter estimates. ",
                   "Check the Mplus output file in '",
                   dirname(file_out), "'.")
        }
    }
    if (!hasEstimates && !any(mDiag$type %in% .notConvergedTypes)) {
        extra[[length(extra) + 1]] <- .diagRow(
            "no_estimates", "error", "Mplus did not produce parameter estimates")
    }

    bounds <- if (hasEstimates && !is.null(model) && !is.null(info)) {
        .activeBounds(model, .mplusLabelEstimates(model, fit), info, data)
    }
    diag <- do.call(.bindDiagnostics, c(list(mDiag), extra, list(bounds), rRows))
    status <- .diagStatus(diag, hasEstimates,
                          llReplicated = .mplusLLReplicated(lines),
                          minCoverage = .mplusMinCoverage(lines))
    list(fit = fit, success = success, warningM = warningM, errorM = errorM,
         diagnostics = diag, status = status)
}

## Short text for the print method
.formatDiagnostics <- function(diag, status, maxRows = 10) {
    yn <- function(x) if (is.na(x)) "NA" else if (isTRUE(x)) "yes" else "no"
    lines <- paste0("Converged: ", yn(status$converged),
                    ", SEs trusted: ", yn(status$seTrusted),
                    ", Admissible: ", yn(status$admissible),
                    ", Best solution replicated: ", yn(status$llReplicated),
                    ", At a bound: ", yn(status$boundary))
    if (!is.na(status$condNumber)) {
        lines <- c(lines, paste0("Condition number: ",
                                 format(signif(status$condNumber, 3))))
    }
    show <- diag[diag$type != "r_warning", , drop = FALSE]
    if (nrow(show) == 0) {
        lines <- c(lines, "No problems detected.")
    } else {
        sev <- c("error", "inadmissible", "local_optimum", "se", "boundary", "info")
        show <- show[order(match(show$severity, sev)), , drop = FALSE]
        txt <- paste0("  [", show$severity, "] ", show$message)
        if (length(txt) > maxRows) {
            txt <- c(txt[seq_len(maxRows)],
                     paste0("  ... and ", length(txt) - maxRows,
                            " more (see $diagnostics)"))
        }
        lines <- c(lines, txt)
    }
    lines
}
