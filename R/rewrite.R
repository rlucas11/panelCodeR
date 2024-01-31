panelcoder <- function(data,
                      panelModel = "starts",
                      program = "lavaan",
                      crossLag = TRUE,
                      stateCors = FALSE,
                      residCors = FALSE,
                      limits = TRUE,
                      stationarity = TRUE,
                      mplusOptions = NULL,
                      mplusOutput = NULL,
                      mplusDirectory = "mplus",
                      lavaanOptions = NULL) {
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
    modelFit <- lavaan::lavaan(model = model, data = data)
    return(list(model,modelFit))
}

                      

