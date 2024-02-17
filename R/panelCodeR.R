.buildModel <- function(data,
                        panelModel = "starts",
                        crossLag = TRUE,
                        stateCors = FALSE,
                        residCors = FALSE,
                        limits = TRUE,
                        stationarity = TRUE,
                        invariance = TRUE
                        ) {
    ## Collect basic info
    info <- getInfo(data)

    if ((panelModel %in% c("starts",
                          "riclpm",
                          "arts",
                          "clpm",
                          "sts")) == FALSE) {
        stop("No model with that name")
    }
    
    if (panelModel == "starts") {
        ar <- TRUE
        trait <- TRUE
        stability <- TRUE
        cl <- TRUE
        state <- TRUE
        traitCors <- TRUE
        arCors <- TRUE
        stateCors <- stateCors
        residCors <- residCors
    }

    if (panelModel == "riclpm") {
        ar <- TRUE
        trait <- TRUE
        stability <- TRUE
        cl <- TRUE
        state <- FALSE
        traitCors <- TRUE
        arCors <- TRUE
        stateCors <- FALSE
        residCors <- residCors
    }

    if (panelModel == "arts") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        cl <- TRUE
        state <- TRUE
        traitCors <- TRUE
        arCors <- TRUE
        stateCors <- stateCors
        residCors <- residCors
    }

    if (panelModel == "clpm") {
        ar <- TRUE
        trait <- FALSE
        stability <- TRUE
        cl <- TRUE
        state <- FALSE
        traitCors <- FALSE
        arCors <- TRUE
        stateCors <- stateCors
        residCors <- residCors
    }

    if (panelModel == "sts") {
        ar <- FALSE
        trait <- TRUE
        stability <- FALSE
        cl <- FALSE
        state <- TRUE
        traitCors <- TRUE
        arCors <- FALSE
        stateCors <- stateCors
        residCors <- residCors
    }

    modelInfo <- .buildTable(info,
        ar = ar,
        trait = trait,
        stability = stability,
        cl = cl,
        state = state,
        traitCors = traitCors,
        arCors = arCors,
        stateCors = stateCors,
        residCors = residCors
    )

    ## Build table and collect parameters
    model <- modelInfo$model

    ## Build final model based on options
    if (ar == TRUE) {
        model <- .constrainStability(model, info)
    }

    ## Constrain cross-lagged paths if necessary
    if (info$gen$yVar == TRUE &
        ar == TRUE &
        crossLag == TRUE) {
        model <- .constrainCl(model, info)
    }

    ## Impose stationarity if requested
    if (stationarity == TRUE &
        ar == TRUE) {
        model <- .arStationarity(model, info)
    }

    ## Constrain state variances
    if (state == TRUE) {
        model <- .constrainStateVar(model, info)
    }

    ## Constrain state correlations
    if (state == TRUE & stateCors == TRUE) {
        model <- .constrainStateCors(model, info)
    }

    if (invariance == TRUE) {
        model <- .constrainLoadings(model, info)
    }

    ## Impose limits on variances and covariances
    ## Eventually change to allow for no correlations
    if (limits == TRUE) {
        model <- .buildLimits(model,
            info,
            ar = ar,
            trait = trait,
            stability = stability,
            cl = cl,
            state = state,
            traitCors = TRUE,
            arCors = TRUE,
            stateCors = stateCors,
            residCors = residCors
        )
    }

    ## Run model in lavaan if requested
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
#'   "riclpm", "clpm", "arts", or "sts". 
#' @param program Program to use to run code. Can be "lavaan" (the default) or
#'   "mplus"
#' @param crossLag Logical value indicating whether to include cross-lagged
#'   paths. Defaults to `TRUE`.
#' @param stateCors Logical value indicating whether to include correlations
#'   between wave-specific state components. Defaults to FALSE.
#' @param residCors Logical value indicating whether to include correlations
#'   between item-specific residual. Defaults to FALSE.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values. Defaults to TRUE.
#' @param stationarity Logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param invariance Logical value indicating whether to constrain loadings for
#'   the same item to be equal across waves. 
#' @param mplusAnalysis Quoted text. Specify ANLYSIS command for mplus. Defaults
#' to "MODEL=NOCOVARIANCES;"
#' @param mplusOutput Quoted text. Specify OUTPUT command for mplus. Defaults to
#'   "stdyx; \\n cinterval; \\n",
#' @param mplusDirectory Quoted text. Specify directory for mplus input and
#'   output files. This directory must already exist before running the command.
#'   Defaults to "mplus".
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator). Not yet implemented
#' @param run Logical value indicating whether to run the model or to just
#'   print and return code.
#' @param ... Additional options passed to lavaan. Default options are
#'   meanstructure=TRUE and missing = 'fiml'.
#' @returns pcObject, which is a list that includes the parameter table used to
#'   create the model, a list of basic information about the model, the model
#'   code (either lavaan or mplus) and the fit object produced by lavaan or
#'   MplusAutomation.
#' @export
panelcoder <- function(data,
                       title = "panelcoder",
                       panelModel = "starts",
                       program = "lavaan",
                       crossLag = TRUE,
                       stateCors = FALSE,
                       residCors = FALSE,
                       limits = TRUE,
                       stationarity = TRUE,
                       invariance = TRUE,
                       mplusAnalysis = NULL,
                       mplusOutput = NULL,
                       mplusDirectory = "mplus",
                       constrainCors = TRUE,
                       run = TRUE,
                       ...
                       ) {
    info <- getInfo(data)
    model <- .buildModel(data,
                         panelModel,
                         crossLag,
                         stateCors,
                         residCors,
                         limits,
                         stationarity,
                         invariance
                         )

    ## Lavaan
    if (program == "lavaan") {
        if (run == TRUE) {
            fit <- lavaan::lavaan(model,
                          data = data,
                          meanstructure = TRUE,
                          missing = 'fiml')
            pcSum <- .summarizeLavaan(fit)
        } else {
            modelCode <- lav2lavaan(model)
            cat(modelCode)
            return(modelCode)
        }
    }
    
    if (program == "mplus") {
        mplusModel <- lav2mplus(model)
        mplusStatement <- MplusAutomation::mplusObject(TITLE = title,
                                      rdata = data,
                                      ANALYSIS = "MODEL=NOCOVARIANCES;",
                                      OUTPUT = "stdyx; cinterval; TECH4 \n",
                                      MODEL = mplusModel)
        if (run == TRUE) {
            fit <- MplusAutomation::mplusModeler(mplusStatement,
                                modelout = paste0(mplusDirectory,
                                                  "/",
                                                  title,
                                                  ".inp"),
                                run = 1)
            pcSum <- .summarizeMplus(info, fit)
        } else {
            cat(mplusModel)
            return(mplusModel)
        }
    }
    ## Get average correlations for each lag
    ## Temporarily only run if using manifest-variable model
    if (info$x$indicators > 1) {
        latent <- 1
    } else {
        latent <- 0
    }

    if (info$gen$yVar == TRUE) {
        if (info$y$indicators > 1) {
            latent <- 1
        }
    }
    
    if (latent == 0) {
        if (info$gen$yVar == TRUE) {
            varNames <- c(info$x$name, info$y$name)
        } else {
            varNames <- info$x$name
        }
        corSummary <- combineCors(data,
                                  info,
                                  program,
                                  fit)
    } else {
        corSummary <- NULL
    }
    pcOutput <- list(pcSum, info, model, fit, corSummary)
    class(pcOutput) <- "pcOutput"
    summary(pcSum)
    return(pcOutput)
}

    