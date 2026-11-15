################################################################################
## Starting values
##
## All forms of the `startVals` argument to panelcoder() resolve to a numeric
## vector with one entry per row of the parameter table (NA = leave the program
## default). .applyStartVals() then writes these into the `ustart` column, but
## only for free parameters that are not dependent on others through an `==`
## constraint. lav2mplus() writes `ustart` of free parameters as `*value`;
## lavaan reads `ustart` directly from the parameter table.
################################################################################

## Base role names. Variable-specific versions add "_x"/"_y" (or "_xy"/"_yx"
## for cross-lagged effects, where "_xy" is the effect of x on y). A base name
## used on its own applies to both variables.
.startValRoleBases <- c(
    "stab", "cl", "ma", "clma",
    "trait", "ar", "innov", "state", "slope",
    "trait_loading", "slope_loading", "loading",
    "resid", "resid_cov",
    "cov_trait", "cov_ar", "cov_innov", "cov_state", "cov_slope",
    "cov_trait_slope", "cov_trait_ar", "cov_slope_ar"
)

.startValRoleNames <- function() {
    c(.startValRoleBases,
      as.vector(outer(.startValRoleBases, c("_x", "_y", "_xy", "_yx"), paste0)))
}

## Labels whose value follows from an equality constraint (lhs of `==`).
.dependentLabels <- function(table) {
    unique(table$lhs[table$op == "=="])
}

## Rows that can receive a starting value
.startValEligible <- function(table) {
    table$free != 0 &
        table$op %in% c("~", "~~", "=~") &
        !(table$label != "" & table$label %in% .dependentLabels(table))
}

## Split a variable name from the parameter table into its parts. Observed
## variables are stem_wave(_indicator); latent variables are
## prefix_stem(_wave).
.parsePcName <- function(name, info) {
    parts <- strsplit(name, "_", fixed = TRUE)[[1]]
    stems <- info$x$name
    if (info$gen$yVar) stems <- c(stems, info$y$name)
    num <- suppressWarnings(as.numeric(parts))
    if (length(parts) >= 2 && !is.na(num[2])) {
        prefix <- "obs"
        stem <- parts[1]
        wave <- num[2]
    } else {
        prefix <- parts[1]
        stem <- if (length(parts) >= 2) parts[2] else NA_character_
        wave <- if (length(parts) >= 3) num[3] else NA_real_
    }
    var <- NA_character_
    if (!is.na(stem) && stem %in% stems) {
        var <- if (stem == info$x$name) "x" else "y"
    }
    list(prefix = prefix, var = var, wave = wave)
}

## Tag each row of the parameter table with a role (NA if none applies)
.startValRoles <- function(table, info) {
    n <- nrow(table)
    roles <- rep(NA_character_, n)
    lhsP <- lapply(table$lhs, .parsePcName, info = info)
    rhsP <- lapply(table$rhs, .parsePcName, info = info)

    ## Is there an autoregressive process (so that later-wave AR variances
    ## are innovations)?
    hasAR <- FALSE
    ## First wave of the AR component for each variable
    firstAR <- c(x = Inf, y = Inf)
    for (r in seq_len(n)) {
        l <- lhsP[[r]]
        rr <- rhsP[[r]]
        if (table$op[r] == "~" && l$prefix == "a" && rr$prefix == "a") {
            hasAR <- TRUE
        }
        if (table$op[r] == "~~" && l$prefix == "i" && !is.na(l$var) &&
            !is.na(l$wave)) {
            firstAR[l$var] <- min(firstAR[l$var], l$wave)
        }
    }

    for (r in seq_len(n)) {
        op <- table$op[r]
        l <- lhsP[[r]]
        rr <- rhsP[[r]]
        role <- NA_character_
        if (op == "~") {
            if (l$prefix == "a" && rr$prefix == "a") {
                role <- if (identical(l$var, rr$var)) paste0("stab_", l$var)
                        else paste0("cl_", rr$var, l$var)
            } else if (l$prefix == "a" && rr$prefix == "i") {
                role <- if (identical(l$var, rr$var)) paste0("ma_", l$var)
                        else paste0("clma_", rr$var, l$var)
            }
        } else if (op == "=~") {
            if (l$prefix == "t" && rr$prefix %in% c("l", "a")) {
                role <- paste0("trait_loading_", l$var)
            } else if (l$prefix == "sl") {
                role <- paste0("slope_loading_", l$var)
            } else if (l$prefix == "l" && rr$prefix == "obs") {
                role <- paste0("loading_", l$var)
            }
        } else if (op == "~~") {
            if (table$lhs[r] == table$rhs[r]) {
                role <- switch(l$prefix,
                    t = paste0("trait_", l$var),
                    i = if (hasAR && !is.na(l$var) && !is.na(l$wave) &&
                            l$wave > firstAR[l$var]) {
                            paste0("innov_", l$var)
                        } else {
                            paste0("ar_", l$var)
                        },
                    s = paste0("state_", l$var),
                    sl = paste0("slope_", l$var),
                    obs = paste0("resid_", l$var),
                    NA_character_)
            } else {
                pre <- sort(c(l$prefix, rr$prefix))
                key <- paste(pre, collapse = "-")
                sameWave <- identical(l$wave, rr$wave)
                role <- switch(key,
                    "obs-obs" = if (identical(l$var, rr$var)) {
                                    paste0("resid_cov_", l$var)
                                } else {
                                    "resid_cov"
                                },
                    "t-t" = "cov_trait",
                    "i-i" = if (sameWave && !identical(l$var, rr$var)) {
                                if (hasAR && l$wave > min(firstAR)) {
                                    "cov_innov"
                                } else {
                                    "cov_ar"
                                }
                            } else {
                                NA_character_
                            },
                    "s-s" = if (sameWave) "cov_state" else NA_character_,
                    "sl-sl" = "cov_slope",
                    "sl-t" = "cov_trait_slope",
                    "i-t" = "cov_trait_ar",
                    "i-sl" = "cov_slope_ar",
                    NA_character_)
            }
        }
        if (!is.na(role) && grepl("_NA", role, fixed = TRUE)) {
            role <- NA_character_
        }
        roles[r] <- role
    }
    roles
}

## Does a user-supplied name match a row role? Generic names (e.g. "stab")
## match both variable-specific roles ("stab_x", "stab_y").
.roleMatches <- function(name, roles) {
    !is.na(roles) &
        (roles == name | sub("_(x|y|xy|yx)$", "", roles) == name)
}

## Does a user-supplied name refer to a role that is more specific?
.roleSpecificity <- function(name) {
    if (grepl("_(x|y|xy|yx)$", name) && !(name %in% .startValRoleBases)) 2 else 1
}

## Is startVals a sequence of start value sets to try in turn?
.isStartValsSeq <- function(startVals) {
    is.list(startVals) &&
        !inherits(startVals, "pcOutput") &&
        is.null(names(startVals))
}

## Short description of a start value set (for the attempt summary)
.describeStartVals <- function(startVals) {
    if (is.null(startVals)) return("default")
    if (is.character(startVals)) return(startVals[1])
    if (inherits(startVals, "pcOutput")) return("previous fit")
    if (is.function(startVals)) return("function")
    "user"
}

## Estimates from a previous panelcoder() result, named by parameter label
.prevFitStartVals <- function(pcOut) {
    fit <- pcOut[[4]]
    oldModel <- pcOut[[3]]
    if (is.null(fit)) {
        stop("The panelcoder() result supplied to startVals has no estimates (the model did not run successfully).",
             call. = FALSE)
    }
    if (inherits(fit, "lavaan")) {
        pt <- lavaan::parTable(fit)
        keep <- pt$label != "" & pt$free > 0
        return(stats::setNames(pt$est[keep], pt$label[keep]))
    }
    est <- fit$results$parameters$unstandardized
    if (is.null(est)) {
        stop("Could not find parameter estimates in the panelcoder() result supplied to startVals.",
             call. = FALSE)
    }
    rows <- which(oldModel$free != 0 & oldModel$label != "" &
                  oldModel$op %in% c("~", "~~", "=~"))
    vals <- rep(NA_real_, length(rows))
    hdr <- toupper(est$paramHeader)
    prm <- toupper(est$param)
    for (k in seq_along(rows)) {
        r <- rows[k]
        lhs <- toupper(gsub(".", "_", oldModel$lhs[r], fixed = TRUE))
        rhs <- toupper(gsub(".", "_", oldModel$rhs[r], fixed = TRUE))
        op <- oldModel$op[r]
        hit <- integer(0)
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
        if (length(hit) > 0) vals[k] <- est$est[hit[1]]
    }
    out <- stats::setNames(vals, oldModel$label[rows])
    out[!is.na(out)]
}

## Latent (wave-level) variance of a variable at one wave, estimated from data.
## With one indicator this is the observed variance. With multiple indicators
## it is the mean covariance of the marker (first) indicator with the others,
## which equals the latent variance when loadings are close to 1.
.waveVariance <- function(data, stem, wave, indicators) {
    if (indicators == 1) {
        col <- paste(stem, wave, sep = "_")
        if (!(col %in% names(data))) return(NA_real_)
        return(stats::var(data[[col]], na.rm = TRUE))
    }
    cols <- paste(stem, wave, seq_len(indicators), sep = "_")
    if (!(cols[1] %in% names(data))) return(NA_real_)
    marker <- data[[cols[1]]]
    covs <- sapply(cols[-1], function(cc) {
        if (!(cc %in% names(data))) return(NA_real_)
        stats::cov(marker, data[[cc]], use = "pairwise.complete.obs")
    })
    psi <- mean(covs, na.rm = TRUE)
    if (!is.finite(psi) || psi <= 0) {
        psi <- .5 * stats::var(marker, na.rm = TRUE)
    }
    psi
}

## Built-in data-based starting values. Each variable's (latent) variance is
## split evenly across the variance components in the model (trait, AR,
## state); stabilities start at .5 (lag 1) or 0 (longer lags); cross-lags,
## moving averages, and covariances start at 0; free loadings at 1.
.dataStartVals <- function(table, info, data, roles, eligible) {
    out <- rep(NA_real_, nrow(table))
    stab <- .5
    vars <- if (info$gen$yVar) c("x", "y") else "x"
    traitOnAR <- any(table$op == "=~" & grepl("^t_", table$lhs) &
                     grepl("^a_", table$rhs))

    for (v in vars) {
        vInfo <- info[[v]]
        stem <- vInfo$name
        Vt <- sapply(vInfo$waves, function(w) {
            .waveVariance(data, stem, w, vInfo$indicators)
        })
        names(Vt) <- vInfo$waves
        Vbar <- mean(Vt, na.rm = TRUE)
        if (!is.finite(Vbar) || Vbar <= 0) next
        Vwave <- function(w) {
            val <- Vt[as.character(w)]
            if (is.na(val) || val <= 0) Vbar else unname(val)
        }

        k <- sum(c(any(roles %in% paste0("trait_", v)),
                   any(roles %in% paste0("ar_", v)),
                   any(roles %in% paste0("state_", v))))
        k <- max(k, 1)

        nElig <- function(role) sum(eligible & roles %in% role)
        poolAR <- nElig(paste0("innov_", v)) == 0 &&
            nElig(paste0("ar_", v)) <= 1
        poolInnov <- nElig(paste0("innov_", v)) <= 1
        poolState <- nElig(paste0("state_", v)) <= 1

        ## Slope scaling: largest fixed slope loading
        slRows <- table$op == "=~" & table$lhs == paste0("sl_", stem) &
            table$free == 0 & !is.na(table$ustart)
        maxLoad <- if (any(slRows)) max(abs(table$ustart[slRows])) else 1
        slWaves <- vInfo$waves

        for (r in which(eligible & !is.na(roles))) {
            role <- roles[r]
            if (!grepl(paste0("_", v, "$"), role)) next
            wave <- .parsePcName(table$lhs[r], info)$wave
            base <- sub(paste0("_", v, "$"), "", role)
            val <- switch(base,
                trait = Vbar / k * (if (traitOnAR) (1 - stab)^2 else 1),
                ar = (if (poolAR) Vbar else Vwave(wave)) / k,
                innov = (if (poolInnov) Vbar else Vwave(wave)) / k *
                    (1 - stab^2),
                state = (if (poolState) Vbar else Vwave(wave)) / k,
                slope = .1 * Vbar / k / max(1, maxLoad^2),
                stab = {
                    lag <- wave - .parsePcName(table$rhs[r], info)$wave
                    if (!is.na(lag) && lag == 1) stab else 0
                },
                ma = 0,
                trait_loading = 1,
                loading = 1,
                slope_loading = {
                    w <- .parsePcName(table$rhs[r], info)$wave
                    (w - min(slWaves)) / max(1, max(slWaves) - min(slWaves))
                },
                resid = {
                    obsVar <- stats::var(data[[table$lhs[r]]], na.rm = TRUE)
                    max(obsVar - Vwave(wave), .1 * obsVar)
                },
                resid_cov = 0,
                NA_real_)
            out[r] <- val
        }
    }

    ## Roles that do not depend on a single variable
    zeroRoles <- c("cl_xy", "cl_yx", "clma_xy", "clma_yx", "resid_cov",
                   "cov_trait", "cov_ar", "cov_innov", "cov_state",
                   "cov_slope", "cov_trait_slope", "cov_trait_ar",
                   "cov_slope_ar")
    out[eligible & roles %in% zeroRoles] <- 0
    out
}

## Resolve a single (non-sequence) startVals specification to a numeric
## vector with one value per row of the table (NA = program default).
.resolveStartVals <- function(startVals, table, info, data) {
    n <- nrow(table)
    out <- rep(NA_real_, n)
    if (is.null(startVals)) return(out)

    eligible <- .startValEligible(table)
    roles <- .startValRoles(table, info)

    if (is.function(startVals)) {
        fArgs <- list(data = data, info = info, table = table)
        fNames <- names(formals(startVals))
        if (!("..." %in% fNames)) fArgs <- fArgs[names(fArgs) %in% fNames]
        res <- do.call(startVals, fArgs)
        if (is.function(res) || .isStartValsSeq(res)) {
            stop("A startVals function must return 'data', a named list or vector of values, or a panelcoder() result.",
                 call. = FALSE)
        }
        return(.resolveStartVals(res, table, info, data))
    }

    if (is.character(startVals)) {
        if (length(startVals) == 1 && startVals == "data") {
            return(.dataStartVals(table, info, data, roles, eligible))
        }
        stop("startVals must be NULL, 'data', a named list or vector of values, a panelcoder() result, a function, or an unnamed list of these.",
             call. = FALSE)
    }

    if (inherits(startVals, "pcOutput")) {
        vals <- .prevFitStartVals(startVals)
        vals <- vals[names(vals) %in% table$label]
        idx <- match(table$label, names(vals))
        hit <- eligible & !is.na(idx) & table$label != ""
        out[hit] <- unname(vals[idx[hit]])
        return(out)
    }

    if (is.numeric(startVals) || is.list(startVals)) {
        nms <- names(startVals)
        if (is.null(nms) || any(nms == "")) {
            stop("Start values supplied by role or label must be named (e.g., list(stab = .5) or c(a_1_2 = .5)).",
                 call. = FALSE)
        }
        vals <- lapply(startVals, function(x) {
            if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
                stop("Each start value must be a single number.", call. = FALSE)
            }
            x
        })
        roleNames <- .startValRoleNames()
        deps <- .dependentLabels(table)
        isRole <- nms %in% roleNames
        unknown <- nms[!isRole & !(nms %in% table$label)]
        if (length(unknown) > 0) {
            warning("Start values ignored for names that are not roles or parameter labels in this model: ",
                    paste(unknown, collapse = ", "), call. = FALSE)
        }
        constrained <- nms[!isRole & nms %in% deps]
        if (length(constrained) > 0) {
            warning("Start values ignored for parameters that are determined by constraints: ",
                    paste(constrained, collapse = ", "), call. = FALSE)
        }
        ## Generic roles first, then specific roles, then exact labels, so
        ## that more specific names override less specific ones
        roleOrder <- nms[isRole]
        roleOrder <- roleOrder[order(sapply(roleOrder, .roleSpecificity))]
        for (nm in roleOrder) {
            hit <- eligible & .roleMatches(nm, roles)
            out[hit] <- vals[[nm]]
        }
        for (nm in nms[!isRole]) {
            hit <- eligible & table$label == nm
            out[hit] <- vals[[nm]]
        }
        return(out)
    }

    stop("startVals must be NULL, 'data', a named list or vector of values, a panelcoder() result, a function, or an unnamed list of these.",
         call. = FALSE)
}

## Write resolved start values into the parameter table
.applyStartVals <- function(table, values) {
    eligible <- .startValEligible(table)
    hit <- eligible & !is.na(values)
    table$ustart[hit] <- signif(values[hit], 4)
    table
}

## Did a panelcoder() attempt converge with usable standard errors?
.startValsConverged <- function(pcOut, program) {
    if (!inherits(pcOut, "pcOutput")) return(FALSE)
    if (is.null(pcOut[[1]]) || !is.null(pcOut[[8]])) return(FALSE)
    fit <- pcOut[[4]]
    if (program == "lavaan") {
        if (!isTRUE(lavaan::lavInspect(fit, "converged"))) return(FALSE)
        pt <- lavaan::parTable(fit)
        if (is.null(pt$se) || !all(is.finite(pt$se[pt$free > 0]))) {
            return(FALSE)
        }
    }
    TRUE
}

## Log-likelihood of a panelcoder() attempt
.startValsLogLik <- function(pcOut, program) {
    if (!inherits(pcOut, "pcOutput") || is.null(pcOut[[4]])) return(NA_real_)
    fit <- pcOut[[4]]
    ll <- tryCatch({
        if (program == "lavaan") {
            unname(lavaan::fitMeasures(fit, "logl"))
        } else {
            fit$results$summaries$LL
        }
    }, error = function(e) NA_real_)
    if (is.null(ll) || length(ll) == 0) NA_real_ else as.numeric(ll[1])
}

## Run panelcoder() once for each start value set in a sequence, returning the
## first converged result (rule = "first") or the converged result with the
## highest log-likelihood (rule = "best").
.runStartValsSeq <- function(args, dots, startValsSeq, rule) {
    if (length(startValsSeq) == 0) {
        stop("startVals is an empty list.", call. = FALSE)
    }
    if (isFALSE(args$run)) {
        stop("Sequences of start values require run = TRUE. Supply a single set of start values to view model code.",
             call. = FALSE)
    }
    program <- args$program
    n <- length(startValsSeq)
    results <- vector("list", n)
    logs <- vector("list", n)
    warns <- character(0)
    attempts <- data.frame(attempt = seq_len(n),
                           startVals = vapply(startValsSeq, .describeStartVals,
                                              character(1)),
                           converged = NA,
                           logLik = NA_real_,
                           error = NA_character_,
                           selected = FALSE,
                           stringsAsFactors = FALSE)
    for (i in seq_len(n)) {
        a <- args
        a["startVals"] <- list(startValsSeq[[i]])
        a$title <- paste0(args$title, "_sv", i)
        res <- NULL
        logs[[i]] <- utils::capture.output(
            res <- withCallingHandlers(
                tryCatch(do.call(panelcoder, c(a, dots)),
                         error = function(e) e),
                warning = function(w) {
                    warns <<- c(warns, conditionMessage(w))
                    invokeRestart("muffleWarning")
                })
        )
        results[[i]] <- res
        if (inherits(res, "error")) {
            attempts$converged[i] <- FALSE
            attempts$error[i] <- conditionMessage(res)
        } else {
            attempts$converged[i] <- .startValsConverged(res, program)
            attempts$logLik[i] <- .startValsLogLik(res, program)
            if (inherits(res, "pcOutput") && !is.null(res[[8]])) {
                attempts$error[i] <- res[[8]]
            } else if (!attempts$converged[i]) {
                attempts$error[i] <- "Did not converge"
            }
        }
        if (rule == "first" && attempts$converged[i]) break
    }

    ok <- which(attempts$converged %in% TRUE)
    if (length(ok) == 0) {
        chosen <- 1
    } else if (rule == "best") {
        ll <- attempts$logLik[ok]
        ## Ties (within .001) go to the earliest attempt
        chosen <- if (all(is.na(ll))) {
            ok[1]
        } else {
            ok[which(ll >= max(ll, na.rm = TRUE) - .001)[1]]
        }
    } else {
        chosen <- ok[1]
    }
    attempts$selected[chosen] <- TRUE

    for (w in unique(warns)) warning(w, call. = FALSE)
    out <- results[[chosen]]
    if (inherits(out, "error")) stop(conditionMessage(out), call. = FALSE)

    cat(logs[[chosen]], sep = "\n")
    if (length(ok) == 0) {
        cat("\nNo start value set produced a converged solution; returning attempt 1.\n")
    } else {
        cat(paste0("\nStart values: used attempt ", chosen, " (",
                   attempts$startVals[chosen], ") of ", n, ".\n"))
    }
    out$startVals <- attempts
    out
}
