# Diagnostics: .mplusDiagnostics() classifies Mplus messages, .activeBounds()
# finds the inequality constraints a solution sits on, and panelcoder() adds
# `diagnostics` and `status` to its result without changing elements 1-8.
#
# Fixtures are trimmed excerpts of real Mplus 8.11 output from models fitted to
# gen_starts() data, except where noted: the first-order-derivative SE warning
# and the THETA warning could not be provoked with simulated data, so those
# excerpts are assembled from the message text in the Mplus 8.11 executable.

write_out <- function(lines) {
    path <- tempfile(fileext = ".out")
    writeLines(lines, path)
    path
}

## --- Fixtures ----------------------------------------------------------------

# Real: STARTS, limits = FALSE, data without state variance
fx_psi <- c(
    "THE MODEL ESTIMATION TERMINATED NORMALLY",
    "",
    "     WARNING:  THE LATENT VARIABLE COVARIANCE MATRIX (PSI) IS NOT POSITIVE",
    "     DEFINITE.  THIS COULD INDICATE A NEGATIVE VARIANCE/RESIDUAL VARIANCE FOR A",
    "     LATENT VARIABLE, A CORRELATION GREATER OR EQUAL TO ONE BETWEEN TWO LATENT",
    "     VARIABLES, OR A LINEAR DEPENDENCY AMONG MORE THAN TWO LATENT VARIABLES.",
    "     CHECK THE TECH4 OUTPUT FOR MORE INFORMATION.",
    "     PROBLEM INVOLVING VARIABLE S_X_1.",
    ""
)

# Assembled from Mplus 8.11 message text
fx_theta <- c(
    "     WARNING:  THE RESIDUAL COVARIANCE MATRIX (THETA) IS NOT POSITIVE DEFINITE.",
    "     THIS COULD INDICATE A NEGATIVE VARIANCE/RESIDUAL VARIANCE FOR AN OBSERVED",
    "     VARIABLE, A CORRELATION GREATER OR EQUAL TO ONE BETWEEN TWO OBSERVED",
    "     VARIABLES, OR A LINEAR DEPENDENCY AMONG MORE THAN TWO OBSERVED VARIABLES.",
    "     CHECK THE RESULTS SECTION FOR MORE INFORMATION.",
    "     PROBLEM INVOLVING VARIABLE X_1_2.",
    ""
)

# Assembled from Mplus 8.11 message text (and the wording check_nonpos()
# matched in real output)
fx_se_untrusted <- c(
    "     THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES MAY NOT BE",
    "     TRUSTWORTHY FOR SOME PARAMETERS DUE TO A NON-POSITIVE DEFINITE",
    "     FIRST-ORDER DERIVATIVE PRODUCT MATRIX.  THIS MAY BE DUE TO THE STARTING",
    "     VALUES BUT MAY ALSO BE AN INDICATION OF MODEL NONIDENTIFICATION.  THE",
    "     CONDITION NUMBER IS       0.647D-12.  PROBLEM INVOLVING THE FOLLOWING PARAMETER:",
    "     Parameter 27, S_X_1",
    ""
)

# Real: factor with two indicators and free loadings
fx_se_failed <- c(
    "THE MODEL ESTIMATION TERMINATED NORMALLY",
    "",
    "     THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES COULD NOT BE",
    "     COMPUTED.  THE MODEL MAY NOT BE IDENTIFIED.  CHECK YOUR MODEL.",
    "     PROBLEM INVOLVING THE FOLLOWING PARAMETER:",
    "     Parameter 6, X_2",
    "",
    "     THE CONDITION NUMBER IS      -0.674D-05.",
    ""
)

# Real: iterations = 3
fx_noconv <- c(
    "RANDOM STARTS RESULTS RANKED FROM THE BEST TO THE WORST FIT FUNCTION VALUES",
    "",
    "Fit function values at local maxima and random start numbers:",
    "",
    "",
    "Unperturbed starting value run did not converge.",
    "",
    "5 starting value run(s) did not converge.",
    "     NO CONVERGENCE.  NUMBER OF ITERATIONS EXCEEDED.",
    ""
)

# Real: STARTS with rstarts = 5
fx_starts <- c(
    "RANDOM STARTS RESULTS RANKED FROM THE BEST TO THE WORST FIT FUNCTION VALUES",
    "",
    "Fit function values at local maxima and random start numbers:",
    "",
    "             21.0167         3",
    "             21.0168         4",
    "             21.0169         1",
    "             21.0171         2",
    "             21.0194         5",
    "",
    "Unperturbed starting value run did not converge.",
    "",
    "",
    "",
    "MODEL FIT INFORMATION"
)

# The same table with the best value found only once
fx_starts_local <- sub("21.0168", "21.0368", fx_starts, fixed = TRUE)
fx_starts_local <- sub("21.0169", "21.0369", fx_starts_local, fixed = TRUE)
fx_starts_local <- sub("21.0171", "21.0371", fx_starts_local, fixed = TRUE)

# The same table with only one converged start
fx_starts_one <- c(fx_starts[1:5],
                   "",
                   "Unperturbed starting value run did not converge.",
                   "",
                   "4 perturbed starting value run(s) did not converge.",
                   "",
                   "MODEL FIT INFORMATION")

# Real: little overlap between x_1 and y_4, COVERAGE = .1
fx_coverage <- c(
    "     PROPORTION OF DATA PRESENT",
    "",
    "",
    "           Covariance Coverage",
    "              X_1           X_2           Y_4",
    "              ________      ________      ________",
    " X_1            0.033",
    " X_2            0.033         1.000",
    " Y_4            0.000         0.013         0.013",
    "",
    "",
    "",
    "     THE COVARIANCE COVERAGE FALLS BELOW THE SPECIFIED LIMIT.",
    "     THE MISSING DATA EM ALGORITHM WILL NOT BE INITIATED.",
    "     CHECK YOUR DATA OR LOWER THE COVARIANCE COVERAGE LIMIT.",
    "",
    "",
    "UNIVARIATE SAMPLE STATISTICS"
)

# Real: rstarts = "4 2"
fx_input_error <- c(
    "*** ERROR in ANALYSIS command",
    "  The STARTS option for TYPE=GENERAL takes only one value.  Final stage",
    "  optimizations are done for all random starts.  Respecify the STARTS",
    "  option with one value.",
    "",
    "",
    "MUTHEN & MUTHEN"
)

## --- Mplus message classification -------------------------------------------

test_that("PSI and THETA warnings are typed with the variable involved", {
    d <- .mplusDiagnostics(fx_psi)
    expect_equal(d$type, "psi_npd")
    expect_equal(d$severity, "inadmissible")
    expect_equal(d$parameter, "S_X_1")
    d <- .mplusDiagnostics(fx_theta)
    expect_equal(d$type, "theta_npd")
    expect_equal(d$parameter, "X_1_2")
})

test_that("untrustworthy SEs are typed, with condition number and parameter", {
    d <- .mplusDiagnostics(fx_se_untrusted)
    expect_equal(d$type, "se_untrusted")
    expect_equal(d$severity, "se")
    expect_equal(d$value, 0.647e-12)
    expect_equal(d$parameter, "Parameter 27, S_X_1")
    expect_true(check_nonpos(write_out(fx_se_untrusted)))
    # Not to be confused with the latent covariance (PSI) warning
    expect_false(check_nonpos(write_out(fx_psi)))
})

test_that("failed SEs pick up the condition number from the next paragraph", {
    d <- .mplusDiagnostics(fx_se_failed)
    expect_equal(d$type, "se_failed")
    expect_equal(d$value, -0.674e-5)
    expect_equal(d$parameter, "Parameter 6, X_2")
    expect_true(check_se(write_out(fx_se_failed)))
})

test_that("nonconvergence reports its variant", {
    d <- .mplusDiagnostics(fx_noconv)
    expect_equal(d$type, "no_convergence")
    expect_match(d$message, "No convergence: number of iterations exceeded")
    expect_true(is.na(.mplusLLReplicated(fx_noconv)))
})

test_that("random-start replication is read from the starts table", {
    expect_true(.mplusLLReplicated(fx_starts))
    expect_equal(nrow(.mplusDiagnostics(fx_starts)), 0)
    expect_false(.mplusLLReplicated(fx_starts_local))
    d <- .mplusDiagnostics(fx_starts_local)
    expect_equal(d$type, "ll_not_replicated")
    expect_equal(d$severity, "local_optimum")
    expect_equal(d$value, 1)
    # No table (no random starts): NA
    expect_true(is.na(.mplusLLReplicated(fx_psi)))
    # Only one start converged: replication cannot be checked, so NA and an
    # info diagnostic rather than ll_not_replicated
    expect_true(is.na(.mplusLLReplicated(fx_starts_one)))
    d <- .mplusDiagnostics(fx_starts_one)
    expect_equal(d$type, "starts_insufficient")
    expect_equal(d$severity, "info")
    expect_match(d$message, "Only 1 of 6 starts converged")
    # Mplus's own message takes precedence
    expect_false(.mplusLLReplicated(c(fx_starts, "",
        "     WARNING:  THE BEST LOGLIKELIHOOD VALUE WAS NOT REPLICATED.  THE",
        "     SOLUTION MAY NOT BE TRUSTWORTHY DUE TO LOCAL MAXIMA.  INCREASE THE")))
})

test_that("coverage warnings carry the minimum coverage", {
    d <- .mplusDiagnostics(fx_coverage)
    expect_equal(d$type, "coverage")
    expect_equal(d$value, 0)
    expect_equal(.mplusMinCoverage(fx_coverage), 0)
    expect_true(check_coverage(write_out(fx_coverage)))
})

test_that("input errors are reported", {
    d <- .mplusDiagnostics(fx_input_error)
    expect_equal(d$type, "mplus_input_error")
    expect_match(d$message, "^\\*\\*\\* ERROR in ANALYSIS command")
})

test_that("unclassified MplusAutomation messages are kept as mplus_warning", {
    fit <- list(results = list(
        warnings = list(c("*** WARNING in MODEL command",
                          "  Variable is uncorrelated with all other variables: X_1")),
        errors = list()))
    d <- .mplusDiagnostics(fx_starts, fit = fit)
    expect_equal(d$type, "mplus_warning")
    expect_match(d$message, "uncorrelated")
    # Messages that are already classified are not repeated
    fit$results$warnings <- list(fx_psi[3:8])
    d <- .mplusDiagnostics(fx_psi, fit = fit)
    expect_equal(d$type, "psi_npd")
})

## --- Problem 1: no overwriting ----------------------------------------------

test_that("every problem is recorded, while [[7]]/[[8]] keep legacy values", {
    lines <- c(fx_se_untrusted, "", fx_coverage)
    d <- .mplusDiagnostics(lines)
    expect_setequal(d$type, c("se_untrusted", "coverage"))
    legacy <- .legacyMplusMessages(d, lines)
    expect_equal(legacy$warningM, "Coverage warning")

    lines <- c(fx_se_failed, "", fx_noconv)
    d <- .mplusDiagnostics(lines)
    expect_setequal(d$type, c("se_failed", "no_convergence"))
    legacy <- .legacyMplusMessages(d, lines)
    expect_equal(legacy$errorM, "No convergence")
    expect_false(legacy$success)

    legacy <- .legacyMplusMessages(.mplusDiagnostics(fx_input_error),
                                   fx_input_error)
    expect_equal(legacy$errorM, check_mplus_error(write_out(fx_input_error)))
})

## --- Problem 2: an R warning does not skip the checks -----------------------

test_that("an R warning during the Mplus run does not skip the checks", {
    skip_if_not_installed("MplusAutomation")
    dir <- withr::local_tempdir()
    lines <- c("INPUT READING TERMINATED NORMALLY", "", fx_se_failed)
    local_mocked_bindings(
        mplusModeler = function(object, modelout, run, ...) {
            out <- sub("\\.inp$", ".out", modelout)
            writeLines(lines, out)
            warning("simulated MplusAutomation warning")
            list(results = list(
                output = lines,
                parameters = list(unstandardized = data.frame(
                    paramHeader = "X_1.WITH", param = "X_2",
                    est = 0, se = 0, est_se = 999, pval = 999))))
        },
        .package = "MplusAutomation"
    )
    d <- make_pc_data(nwaves = 4)
    out <- NULL
    utils::capture.output(
        out <- panelcoder(d, panelModel = "clpm", program = "mplus",
                          mplusDirectory = dir, title = "mocked")
    )
    expect_s3_class(out, "pcOutput")
    expect_equal(out[[8]], "Standard errors could not be computed")
    expect_equal(out[[7]], "simulated MplusAutomation warning")
    expect_true(all(c("r_warning", "se_failed") %in% out$diagnostics$type))
    expect_false(out$status$seTrusted)
})

## --- Active bounds -----------------------------------------------------------

test_that("variance and correlation bounds are detected from estimates", {
    d <- make_pc_data(nwaves = 5, bivariate = TRUE)
    info <- getInfo(d)
    model <- .buildModel(data = d, panelModel = "starts")
    labs <- model$label[model$free != 0 & model$label != "" &
                        model$op %in% c("~", "~~", "=~")]
    le <- data.frame(label = labs, est = .5, se = .1)
    le$est[le$label %in% c("x_tVar", "y_tVar")] <- 1
    le$est[le$label == "sx1"] <- 0
    le$est[le$label == "cov_txty"] <- .99
    b <- .activeBounds(model, le, info, d)
    expect_setequal(b$parameter, c("sx1", "cov_txty"))
    expect_true("state variance of x at 0" %in% b$message)
    expect_true("stable-trait correlation at the .99 bound" %in% b$message)
    expect_equal(b$value[b$parameter == "cov_txty"], .99)

    # Negative bound, and nothing flagged in the interior
    le$est[le$label == "cov_txty"] <- -.995
    le$est[le$label == "sx1"] <- .2
    b <- .activeBounds(model, le, info, d)
    expect_equal(b$parameter, "cov_txty")
    expect_match(b$message, "at the -.99 bound")
    le$est[le$label == "cov_txty"] <- .5
    expect_equal(nrow(.activeBounds(model, le, info, d)), 0)
})

test_that("status flags summarize the diagnostics", {
    d <- .bindDiagnostics(.mplusDiagnostics(fx_se_untrusted),
                          .diagRow("boundary", "boundary", "state variance of x at 0"))
    s <- .diagStatus(d, hasEstimates = TRUE, llReplicated = TRUE)
    expect_true(s$converged)
    expect_false(s$seTrusted)
    expect_true(s$admissible)
    expect_true(s$boundary)
    expect_equal(s$condNumber, 0.647e-12)
    s <- .diagStatus(.mplusDiagnostics(fx_noconv), hasEstimates = TRUE)
    expect_false(s$converged)
    s <- .diagStatus(.emptyDiagnostics(), hasEstimates = FALSE)
    expect_false(s$converged)
    expect_true(is.na(s$seTrusted))
})

## --- End to end (lavaan) -----------------------------------------------------

no_state_data <- function() {
    set.seed(1)
    d <- gen_starts(n = 500, nwaves = 5, xr = 0, yr = 0)
    names(d) <- paste(rep(c("x", "y"), each = 5), rep(1:5, 2), sep = "_")
    d
}

test_that("a STARTS fit without state variance reports the boundary (lavaan)", {
    skip_on_cran()
    skip_if_not_installed("lavaan")
    out <- NULL
    utils::capture.output(out <- panelcoder(no_state_data(), panelModel = "starts"))
    expect_null(out[[7]])
    expect_null(out[[8]])
    expect_true(out$status$converged)
    expect_false(out$status$seTrusted)
    expect_true(out$status$boundary)
    expect_true("state variance of x at 0" %in%
                out$diagnostics$message[out$diagnostics$type == "boundary"])
    res <- .collectResults(out)
    expect_match(res$boundary, "state variance of x at 0")
    expect_false(res$seTrusted)
    expect_match(res$diagnostics, "boundary")
    expect_true(all(c("warning", "error", "trait.cor", "ar.cor", "state.cor")
                    %in% names(res)))
    # The summary prints a diagnostics section
    printed <- utils::capture.output(print(out[[1]]))
    expect_true(any(grepl("^Diagnostics", printed)))
    expect_true(any(grepl("state variance of x at 0", printed)))
})

test_that("all lavaan warnings are recorded, not only the last", {
    skip_if_not_installed("lavaan")
    rw <- c("first warning", "second warning")
    d <- .lavaanDiagnostics(NULL, rw)
    expect_equal(d$message, rw)
    expect_equal(unique(d$type), "r_warning")
})

test_that(".collectResults() keeps its old columns for failed runs", {
    res <- .collectResults(list(NULL))
    expect_true(all(c("trait", "ar", "st", "stab", "aic", "bic", "warning",
                      "error", "converged", "seTrusted", "admissible",
                      "llReplicated", "boundary", "condNumber",
                      "diagnostics") %in% names(res)))
    expect_equal(nrow(res), 1)
    expect_true(is.na(res$converged))
})

## --- End to end (Mplus) ------------------------------------------------------

test_that("a STARTS fit without state variance reports the boundary (Mplus)", {
    skip_on_cran()
    skip_if_not_installed("MplusAutomation")
    skip_if(Sys.which("mplus") == "", "Mplus is not installed")
    dir <- withr::local_tempdir()
    out <- NULL
    utils::capture.output(out <- panelcoder(no_state_data(), panelModel = "starts",
                                            program = "mplus",
                                            mplusDirectory = dir,
                                            title = "nostate"))
    expect_null(out[[8]])
    expect_true(out$status$converged)
    expect_false(out$status$seTrusted)
    expect_true(out$status$llReplicated)
    expect_true("state variance of x at 0" %in% out$diagnostics$message)
    expect_true("se_undefined" %in% out$diagnostics$type)
})

## --- SEs printed as 0.000 or as asterisks ------------------------------------

# Label-estimate table as .mplusLabelEstimates() returns it
se_table <- function(se, stdse = NULL, est = .5) {
    le <- data.frame(label = paste0("p", seq_along(se)), lhs = "a", op = "~",
                     rhs = "b", est = est, se = se, stringsAsFactors = FALSE)
    if (!is.null(stdse)) le$stdse <- stdse
    le
}

test_that("Mplus SEs of 0 with a nonzero standardized SE are rounding", {
    d <- .undefinedSEs(se_table(se = 0, stdse = .019, est = .003), NULL,
                       rounded = TRUE)
    expect_equal(d$type, "se_rounded")
    expect_equal(d$severity, "info")
    expect_equal(d$parameter, "p1")
    expect_true(.diagStatus(d, hasEstimates = TRUE)$seTrusted)
})

test_that("undefined Mplus SEs are still flagged", {
    # 999, 0 with standardized SE 0, 0 with standardized SE missing, and NA
    le <- se_table(se = c(999, 0, 0, NA, .1), stdse = c(.1, 0, NA, .1, .1))
    d <- .undefinedSEs(le, NULL, rounded = TRUE)
    expect_equal(d$type, rep("se_undefined", 4))
    expect_equal(d$parameter, c("p1", "p2", "p3", "p4"))
    expect_false(.diagStatus(d, hasEstimates = TRUE)$seTrusted)
})

test_that("without standardized results, small estimates with SE 0 are skipped", {
    d <- .undefinedSEs(se_table(se = c(0, 0), est = c(.003, .5)), NULL,
                       rounded = TRUE)
    expect_equal(d$type, c("se_undefined", "se_rounded"))
    expect_equal(d$parameter, c("p2", "p1"))
})

test_that("lavaan SEs of 0 are always undefined", {
    d <- .undefinedSEs(se_table(se = 0, est = .003), NULL)
    expect_equal(d$type, "se_undefined")
})

test_that("asterisk SEs do not break the summary and are flagged (Mplus)", {
    skip_if_not_installed("MplusAutomation")
    fixture <- test_path("fixtures", "mplus-asterisk-se.out")
    dir <- withr::local_tempdir()
    local_mocked_bindings(
        mplusModeler = function(object, modelout, run, ...) {
            out <- sub("\\.inp$", ".out", modelout)
            file.copy(fixture, out, overwrite = TRUE)
            list(results = MplusAutomation::readModels(out, quiet = TRUE))
        },
        .package = "MplusAutomation"
    )
    d <- make_pc_data(nwaves = 3, nobs = 300)
    out <- NULL
    utils::capture.output(
        out <- panelcoder(d, panelModel = "clpm", program = "mplus",
                          mplusDirectory = dir, title = "asterisk")
    )
    expect_s3_class(out, "pcOutput")
    expect_false(is.null(out[[1]]))
    expect_null(out[[8]])
    se <- out[[4]]$results$parameters$unstandardized$se
    expect_type(se, "double")
    expect_equal(sum(is.na(se)), 1)
    und <- out$diagnostics[out$diagnostics$type == "se_undefined", ]
    expect_equal(nrow(und), 1)
    expect_match(und$message, "I_X_1|i_x_1", ignore.case = TRUE)
    expect_match(und$message, "missing")
    expect_true(out$status$converged)
    expect_false(out$status$seTrusted)
})

# Bivariate data with Y on a scale 100 times larger, so the X-on-Y cross-lag
# has an unstandardized SE that prints as 0.000
rescaled_data <- function() {
    set.seed(2)
    d <- gen_starts(n = 500, nwaves = 5, xr = .5, yr = .5)
    names(d) <- paste(rep(c("x", "y"), each = 5), rep(1:5, 2), sep = "_")
    d[, 6:10] <- d[, 6:10] * 100
    d
}

test_that("SEs that round to 0.000 do not make SEs untrusted (Mplus)", {
    skip_on_cran()
    skip_if_not_installed("MplusAutomation")
    skip_if(Sys.which("mplus") == "", "Mplus is not installed")
    dir <- withr::local_tempdir()
    out <- NULL
    utils::capture.output(out <- panelcoder(rescaled_data(), panelModel = "clpm",
                                            program = "mplus",
                                            mplusDirectory = dir,
                                            title = "rounded"))
    expect_true(out$status$converged)
    expect_true(out$status$seTrusted)
    expect_false("se_undefined" %in% out$diagnostics$type)
    expect_true("se_rounded" %in% out$diagnostics$type)

    # A start value sequence stops at the first attempt
    utils::capture.output(out <- panelcoder(rescaled_data(), panelModel = "clpm",
                                            program = "mplus",
                                            mplusDirectory = dir,
                                            title = "roundedsv",
                                            startVals = list(NULL, "data")))
    expect_true(out$startVals$selected[1])
    expect_true(out$startVals$seTrusted[1])
    expect_true(is.na(out$startVals$converged[2]))
})
