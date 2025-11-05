.buildModel <- function(data,
                        panelModel = "starts",
                        predetermined = FALSE,
                        crossLag = TRUE,
                        lags = 1,
                        traitCors = TRUE,
                        arCors = TRUE,
                        stateCors = FALSE,
                        residCors = FALSE,
                        limits = TRUE,
                        stationarity = "paths",
                        constrainState = TRUE,
                        invariance = TRUE,
                        residVar = FALSE,
                        ma = FALSE,
                        clma = FALSE,
                        slope = "none",
                        state = FALSE,
                        measurement = FALSE
                        ) {
    ## Collect basic info
    info <- getInfo(data)

    ## Stop if problematic specification
    if (info$x$indicators == 1 & panelModel == "measurement") {
        stop("Measurement model not run when there is only 1 indicator")
    }

    if ((panelModel %in% c("starts",
                          "riclpm", "start",
                          "arts",
                          "clpm", "art",
                          "sts",
                          "dpm_c", "dpm_p",
                          "gclm",
                          "lgcm",
                          "alt", "lcmsr",
                          "measurement")) == FALSE) {
        stop("No model with that name")
    }

    if ((slope %in% c("linear",
                      "basis",
                      "centered",
                      "none"
                      )) == FALSE) {
        stop("Use 'linear', 'basis', or 'centered' for slope option, or specify 'none'.")
    }

    if (panelModel == "lgcm" & slope == "none") {
        stop("Need to specify 'linear', 'centered', or 'basis' for slope option when using LGCM")
    }
    

    if (panelModel == "dpm_c" & slope != "none") {
        stop("Bivariate Constrained DPM model is not yet correctly specified when slopes are included. Either use the Predetermined DPM or manually modify code to get desired results.")
    }

    if (panelModel == "starts") {
        ar <- TRUE
        trait <- TRUE
        stability <- TRUE
        crossLag <- crossLag
        state <- TRUE
        traitCors <- traitCors
        arCors <- arCors
        stateCors <- stateCors
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- slope
    }

    if (panelModel == "riclpm" | panelModel == "start") {
        ar <- TRUE
        trait <- TRUE
        stability <- TRUE
        crossLag <- crossLag
        state <- FALSE
        traitCors <- traitCors
        arCors <- arCors
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- slope
    }

    if (panelModel == "lcmsr") {
        ar <- TRUE
        trait <- TRUE
        stability <- TRUE
        crossLag <- crossLag
        state <- state
        traitCors <- traitCors
        arCors <- arCors
        stateCors <- stateCors
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- "linear"
    }

    if (panelModel == "arts") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        crossLag <- crossLag
        state <- TRUE
        traitCors <- FALSE
        arCors <- arCors
        stateCors <- stateCors
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- slope
    }

    if (panelModel == "clpm" | panelModel == "art") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        crossLag <- crossLag
        state <- FALSE
        traitCors <- FALSE
        arCors <- arCors
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- slope
    }
    

    if (panelModel == "sts") {
        ar <- FALSE
        trait <- TRUE
        stability <- FALSE
        crossLag <- FALSE
        state <- TRUE
        traitCors <- traitCors
        arCors <- FALSE
        stateCors <- stateCors
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- FALSE
        clma <- FALSE
        slope <- slope
    }

    if (panelModel == "dpm_c") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        crossLag <- TRUE
        state <- state
        traitCors <- FALSE
        arCors <- arCors
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- TRUE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- slope
    }

    if (panelModel == "dpm_p") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        crossLag <- TRUE
        state <- state
        traitCors <- FALSE
        arCors <- arCors
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- TRUE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- slope
    }

    if (panelModel == "alt") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        crossLag <- TRUE
        state <- state
        traitCors <- FALSE
        arCors <- arCors
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- TRUE
        gclm <- FALSE
        ma <- ma
        clma <- clma
        slope <- "linear"
    }
        
    if (panelModel == "gclm") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        crossLag <- TRUE
        state <- state
        traitCors <- TRUE
        arCors <- arCors
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- TRUE
        ma <- ma
        clma <- clma
        slope <- slope
    }

    if (panelModel == "lgcm") {
        ar <- FALSE
        trait <- TRUE
        stability <- FALSE
        crossLag <- FALSE
        state <- TRUE
        traitCors <- TRUE
        arCors <- FALSE
        stateCors <- stateCors
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        stationarity <- "none"
        ma <- FALSE
        clma <- FALSE
        slope <- slope
    }

    if (panelModel == "measurement") {
        ar <- TRUE
        trait <- FALSE
        stability <- FALSE
        crossLag <- FALSE
        state <- FALSE
        traitCors <- FALSE
        arCors <- FALSE
        stateCors <- FALSE
        residCors <- residCors
        dpm_c <- FALSE
        dpm_p <- FALSE
        gclm <- FALSE
        ma <- FALSE
        clma <- FALSE
        slope <- "none"
        measurement <- TRUE
        stationarity <- "none"
    }


    modelInfo <- .buildTable(info,
        ar = ar,
        trait = trait,
        state = state,
        stability = stability,
        lags = lags,
        crossLag = crossLag,
        traitCors = traitCors,
        arCors = arCors,
        stateCors = stateCors,
        residCors = residCors,
        predetermined = predetermined,
        dpm_c = dpm_c,
        dpm_p = dpm_p,
        gclm = gclm,
        ma = ma,
        clma = clma,
        slope = slope,
        measurement = measurement
    )

    ## Build table and collect parameters
    model <- modelInfo$model

    ## Build final model based on options
    if (ar == TRUE & (stationarity ==  "paths" | stationarity == "full")) {
        model <- .constrainStability(model, info, lags)
    }

    if (ma == TRUE) {
        model <- .constrainMa(model, info)
    }

    if (clma == TRUE) {
        model <- .constrainClMa(model, info)
    }



    ## Constrain cross-lagged paths if necessary
    if (info$gen$yVar == TRUE &
        ar == TRUE &
        measurement == FALSE) {
        if (stationarity == "paths" | stationarity == "full") {
            if (crossLag == FALSE) {
                zero <- TRUE
            } else {
                zero <- FALSE
            }
            model <- .constrainCl(model, info, zero, lags)
        } else if (stationarity == "none" & crossLag == FALSE) {
            zero <- TRUE
            model <- .constrainCl(model, info, zero, lags)
        }
    }

    ## Impose stationarity if requested
    if (stationarity == "full" &
        ar == TRUE) {
        if (arCors == FALSE) {
            zero <- TRUE
        } else {
            zero <- FALSE
        }
        model <- .arStationarity(model, info, zero=zero)
    }

    ## Constrain state variances
    if (state == TRUE & constrainState == TRUE) {
        model <- .constrainStateVar(model, info)
    }

    ## Constrain state correlations
    if (state == TRUE & stateCors == TRUE) {
        model <- .constrainStateCors(model, info)
    }

    if (invariance == TRUE) {
        model <- .constrainLoadings(model, info)
    }

    if (residVar == TRUE) {
        model <- .constrainResidVar(model, info)
    }

    if (dpm_c == TRUE) {
        model <- .constrainDpmLoadings(model, info)
    }

    ## Impose limits on variances and covariances
    ## Eventually change to allow for no correlations
    if (limits == TRUE) {
        model <- .buildLimits(parTable = model,
                              info,
                              panelModel = panelModel,
                              ar = ar,
                              trait = trait,
                              stability = stability,
                              crossLag = crossLag,
                              state = state,
                              constrainState = constrainState,
                              stationarity = stationarity,
                              traitCors = TRUE,
                              arCors = TRUE,
                              stateCors = stateCors,
                              residCors = residCors
                              )
    }
    

    return(model)
}

################################################################################
## panelcoder command
################################################################################

#' Build and Run Lavaan and Mplus Code for Various Panel Models
#'
#' `panelcoder()` produces lavaan and mplus code for variations of the general
#' Stable Trait, Autoregressive Trait, State Model. Univariate and bivariate
#' versions can be specified, and various reduced versions of the model can be
#' run. Other options (such as whether to include lagged paths) can also be
#' selected.
#'
#' @param data Dataframe with appropriately named variables. Names should
#'   include the name stem and wave, separated by "_". If there are multiple
#'   indicators, names should have an additional number (identifying indicator
#'   number), again separated with "_". The number of indicators must be the
#'   same for different waves, but different variables can have different
#'   numbers of indicators. Mplus requires all variable names to be less than
#'   8 characters, so variable names (including the wave and indicator indexes)
#'   should be very short. A warning is provided if these exceed 8 characters.
#' @param title Title of analysis for mplus
#' @param panelModel Specific model to run. Can be "starts"(the default),
#'   "riclpm", "clpm", "arts", "sts", "dpm", "gclm", "alt", "lcmsr", or
#'   "measurement".
#' @param predetermined Logical value indicating whether to use a
#'   "predetermined" version of a STARTS variant (see Andersen, 2022). This
#'   option is not used for the DPM or GCLM. 
#' @param program Program to use to run code. Can be "lavaan" (the default) or
#'   "mplus"
#' @param crossLag Logical value indicating whether to include cross-lagged
#'   paths. Defaults to `TRUE`.
#' @param lags Numeric value indicating the number of lags to include for
#'   stability coefficients and cross-lagged paths. Defaults to 1. Note that
#'   lags greater than 1 cannot yet be included with random intercepts. Also,
#'   full stationarity is not yet implemented with lags greater than 1.
#' @param ma Logical value indicating whether to include moving average
#'   components in the GCLM. Defaults to `FALSE`.
#' @param clma Logical value indicating whether to include cross-lagged moving
#'   average components in the GCLM. Defaults to `FALSE`.
#' @param traitCors Logical value indicating whether to include correlations
#'   between stable-trait/random-intercept components. Defaults to TRUE.
#' @param arCors Logical value indicating whether to include correlations
#'   between wave-specific AR components. Defaults to TRUE.
#' @param stateCors Logical value indicating whether to include correlations
#'   between wave-specific state components. Defaults to FALSE.
#' @param residCors Logical value indicating whether to include correlations
#'   between item-specific residual. Defaults to FALSE.
#' @param residVar Logical value indicating whether to constrain residual
#'   variance to be equal across waves when there are multiple indicators.
#' @param slope String variable indicating what type of slope to specify in
#'   models that include a slope. Can be "none," (the default)"linear, ""basis,"
#'   or "centered."
#' @param state Logical value indicating whether to include a state component.
#'   This is automatically included in some models (e.g., STARTS), but this
#'   option allows for the inclusion of a state component in additional models
#'   including the DPM and GCLM.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values. Defaults to TRUE.
#' @param rstarts Numeric value indicating the number of random starts in
#'   mplus. This is often necessary to avoid local maximum when limiting variance
#'   to be positive, thus the default is to set it to 5 when limits are also set.
#'   If you wish to remove this option, set `rstarts` to NULL.
#' @param stationarity Logical value indicating whether and how to impose
#'   stationarity in the autoregressive process. Defaults to "paths" which
#'   constrains the stability and cross-lagged paths to be equal. Can also be
#'   set to "full" which also constrains variance in the autoregressive
#'   process to be equal across waves or "none" which removes all constraints.
#' @param constrainState Logical value indicating whether to constrain state
#'   variances to be equal across waves. Defaults to TRUE.
#' @param invariance Logical value indicating whether to constrain loadings for
#'   the same item to be equal across waves. 
#' @param mplusAnalysis Quoted text. Specify ANLYSIS command for mplus. Defaults
#'   to "MODEL=NOCOVARIANCES;\\nCOVERAGE=.001;". If you change this, including
#'   "MODEL=NOCOVARIANCES" is highly recommended given the way the model is
#'   specified. 
#' @param mplusOutput Quoted text. Specify OUTPUT command for mplus. Defaults to
#'   "stdyx; \\n cinterval; \\n",
#' @param mplusDirectory Quoted text. Specify directory for mplus input and
#'   output files. This directory must already exist before running the command.
#'   Defaults to "mplus".
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator). Not yet implemented
#' @param run Logical value indicating whether to run the model or to just
#'   print and return code. Defaults to `TRUE`.
#' @param ... Additional options passed to lavaan. Default options are
#'   meanstructure=TRUE, missing = 'fiml', int.ov.free=TRUE, and
#'   int.lv.free=FALSE.
#' @returns pcObject, which is a list that includes the parameter table used to
#'   create the model, a list of basic information about the model, the actual
#'   and implied stability coefficients used for plotting, the model code and
#'   the fit object produced by lavaan or MplusAutomation.
#' @export
panelcoder <- function(data,
                       title = "panelcoder",
                       panelModel = "starts",
                       predetermined = FALSE,
                       program = "lavaan",
                       crossLag = TRUE,
                       lags = 1,
                       ma = FALSE,
                       clma = FALSE,
                       traitCors = TRUE,
                       arCors = TRUE,
                       stateCors = FALSE,
                       residCors = FALSE,
                       residVar = FALSE,
                       slope = "none",
                       state = FALSE,
                       limits = TRUE,
                       rstarts = 5,
                       stationarity = "paths",
                       constrainState = TRUE,
                       invariance = TRUE,
                       mplusAnalysis = NULL,
                       mplusOutput = NULL,
                       mplusDirectory = "mplus",
                       constrainCors = TRUE,
                       run = TRUE,
                       ...
                       ) {
    info <- getInfo(data)

    ## Set impossible options
    if (panelModel == "clpm") {
        stateCors <- FALSE
        traitCors <- FALSE
    }

    if (panelModel == "riclpm") {
        stateCors <- FALSE
    }

    if (panelModel == "arts") {
        traitCors <- FALSE
    }

    if (panelModel == "sts") {
        arCors <- FALSE
    }

    ## Set other conflicting options
    if (limits == FALSE) {
        rstarts == NULL
    }

    ## Check if univariate gclm with clma
    if (panelModel == "gclm" & info$gen$yVar == FALSE & clma == TRUE) {
        clma <- FALSE
        warning("Cross-lagged moving averages not possible for univariate models. CLMA set to FALSE")
    }
    
    
    ## Check for phantom variables
    if (length(info$x$waves) != length(info$x$actualWaves) &
        stationarity != "full") {
        stop("Can't have phantom variables when full stationarity is not set",
             call. = FALSE)
    }

    if (info$gen$yVar == TRUE) {
        if (length(info$y$waves) != length(info$y$actualWaves) &
            stationarity != "full") {
            stop("Can't have phantom variables when full stationarity is not set",
                 call. = FALSE)
        }
    }

    if (panelModel == "dpm" | panelModel == "gclm") {
        if (length(info$x$waves) != length(info$x$actualWaves)) {
            stop("Can't have phantom variables when fitting the dynamic panel model or generalized cross-lagged panel model",
                 call. = FALSE)
        }
    }
    if (panelModel == "dpm_c" |panelModel == "dpm_p" | panelModel == "gclpm") {
        if (info$gen$yVar == TRUE) {    
            if (length(info$y$waves) != length(info$y$actualWaves)) {
                stop("Can't have phantom variables when fitting the dynamic panel model or generalized cross-lagged panel model",
                     call. = FALSE)
            }
        }
    }

    if ((panelModel == "alt" | panelModel == "lcmsr") & slope != "linear") {
        slope = "linear"
        warning("Slope set to linear")
    }

    if (lags > 1 &
        panelModel %in% c("starts", "riclpm", "sts", "gclm", "lgcm", "alt", "lcmsr",
                          "dpm_c", "dpm_p")) {
        stop("Lags greater than 1 not implemented when there is a random intercept. Choose a different model.")
    }

    if (lags > 1 & stationarity == "full") {
        warning("Stationarity not implemented with lags greater than 1. Stationarity set to paths.")
        stationarity = "paths"
    }

    
    model <- .buildModel(data = data,
                         panelModel = panelModel,
                         predetermined = predetermined,
                         crossLag = crossLag,
                         lags = lags,
                         ma = ma,
                         clma = clma,
                         stateCors = stateCors,
                         residCors = residCors,
                         arCors = arCors,
                         slope = slope,
                         state = state,
                         limits = limits,
                         stationarity = stationarity,
                         constrainState = constrainState,
                         invariance = invariance,
                         residVar = residVar
                         )

    ## Lavaan
    if (program == "lavaan") {
        modelCode <- lav2lavaan(model)
        if (run == TRUE) {
            fit <- NULL
            warningM <- NULL
            errorM <- NULL
            output <- tryCatch(
                {
                    fit <- lavaan::lavaan(model,
                        data = data,
                        meanstructure = TRUE,
                        missing = "fiml",
                        int.ov.free = TRUE,
                        int.lv.free = FALSE,
                        ...
                        )
                    list(
                        fit = fit,
                        success = TRUE,
                        warningM = NULL,
                        errorM = NULL
                    )
                },
                warning = function(w) {
                    if (is.null(fit)) {
                        success <- FALSE
                    }
                    list(
                        fit = fit,
                        success = success,
                        warningM = conditionMessage(w),
                        errorM = errorM
                    )
                },
                error = function(e) {
                    list(
                        fit = fit,
                        success = FALSE,
                        warningM = warningM,
                        errorM = conditionMessage(e)
                    )
                }
            )
            if (output$success == TRUE) {
                pcSum <- .summarizeLavaan(panelModel,
                    info = info,
                    fitObject = output$fit,
                    crossLag = crossLag,
                    ma = ma,
                    clma = clma,
                    traitCors = traitCors,
                    arCors = arCors,
                    stateCors = stateCors,
                    residCors = residCors,
                    slope = slope,
                    stationarity = stationarity
                )
            } else {
                pcOutput <- list(
                    NULL,
                    info,
                    model,
                    NULL,
                    NULL,
                    modelCode,
                    output$warningM,
                    output$errorM
                )
                if (!is.null(output$errorM)) {
                    print(paste0("lavaan error: ", output$errorM))
                }
                if (!is.null(output$warningM)) {
                    print(paste0("lavaan warning: ", output$warningM))
                }
                class(pcOutput) <- "pcOutput"
                return(pcOutput)
            }
        } else {
            cat(modelCode)
            return(modelCode)
        }
    }
    
    if (program == "mplus") {
        if (is.null(mplusAnalysis)) {
            mplusAnalysis <- "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;"
        }
        if (!is.null(rstarts)) {
            mplusAnalysis <- paste0(mplusAnalysis,
                                    "STARTS=",
                                    rstarts,
                                    ";\n")
        }
        if (is.null(mplusOutput)) {
            mplusOutput <- "stdyx; cinterval; TECH4; \n"
        }
        
        modelCode <- lav2mplus(model)
        mplusStatement <- MplusAutomation::mplusObject(TITLE = title,
                                      rdata = data,
                                      ANALYSIS = mplusAnalysis,
                                      OUTPUT = mplusOutput,
                                      MODEL = modelCode)
        if (run == TRUE) {
            fit <- NULL
            warningM <- NULL
            errorM <- NULL
            success <- TRUE
            output <- tryCatch(
            {
                file_stem <- paste0(
                    mplusDirectory,
                    "/",
                    title
                )
                file_inp <- paste0(
                    file_stem,
                    ".inp"
                )
                file_out <- paste0(
                    file_stem,
                    ".out"
                )
                fit <- MplusAutomation::mplusModeler(mplusStatement,
                    modelout = file_inp,
                    run = 1
                )
                nonpos <- check_nonpos(file_out)
                if (nonpos) {
                    warningM <- "Non-positive definite matrix"
                }
                noSe <- check_se(file_out)
                if (noSe) {
                    success <- FALSE
                    errorM <- "Standard errors could not be computed"
                }
                noConv <- check_convergence(file_out)
                if (noConv) {
                    success <- FALSE
                    errorM <- "No convergence"
                }
                list(
                    fit = fit,
                    success = success,
                    warningM = warningM,
                    errorM = errorM
                    )
                },
                warning = function(w) {
                    warningM = conditionMessage(w)
                    list(
                        fit = fit,
                        success = TRUE,
                        warningM = warningM,
                        errorM = errorM
                    )
                },
            error = function(e) {
                conditionMessage(e)
                list(
                    fit = fit,
                    success = FALSE,
                    warningM = warningM,
                    errorM = errorM
                )
            }
            )

            if (output$success == TRUE) {
                pcSum <- .summarizeMplus(panelModel,
                    info,
                    output$fit,
                    crossLag = crossLag,
                    ma = ma,
                    clma = clma,
                    traitCors = traitCors,
                    arCors = arCors,
                    stateCors = stateCors,
                    residCors = residCors,
                    slope = slope,
                    stationarity = stationarity
                )
            } else {
                pcOutput <- list(
                    NULL,
                    info,
                    model,
                    NULL,
                    NULL,
                    modelCode,
                    output$warningM,
                    output$errorM
                )
                if (!is.null(output$errorM)) {
                    print(paste0("Mplus error: ", output$errorM))
                }
                if (!is.null(output$warningM)) {
                    print(paste0("Mplus warning: ", output$warningM))
                }
                class(pcOutput) <- "pcOutput"
                return(pcOutput)
            }
        } else {
            cat(modelCode)
            return(modelCode)
        }
    }
    ## Get average correlations for each lag
    ## Temporarily only run if using manifest-variable model
    if (info$x$indicators > 1) {
        latent <- TRUE
    } else {
        latent <- FALSE
    }

    if (info$gen$yVar == TRUE) {
        if (info$y$indicators > 1) {
            latent <- TRUE
        }
    }
    
    if (info$gen$yVar == TRUE) {
        varNames <- c(info$x$name, info$y$name)
    } else {
        varNames <- info$x$name
    }
    corSummary <- combineCors(data,
                              info,
                              program,
                              output$fit,
                              latent)
    
    pcOutput <- list(pcSum, info, model, fit, corSummary, modelCode, output$warningM, output$errorM)
    class(pcOutput) <- "pcOutput"
    print(pcSum)
    if (!is.null(output$warningM)) {
        if (program == "lavaan") {
            print(paste0("lavaan warning: ", output$warningM))
        }
        if (program == "mplus") {
            print(paste0("Mplus warning: ", output$warningM))
        }
    }
    return(pcOutput)
}


