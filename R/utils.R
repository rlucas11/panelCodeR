#' Get info about dataframe
#' @param df dataframe
#' @noRd
getInfo <- function(df) {
    ## Function to get information about names, waves, and indicators from dataframe
    ## Returns errors if data are not structured and named correctly
    ## Warns if length of names exceeds 8 characters (problem for Mplus)
    dfNames <- names(df)

    ## Check length of names for Mplus
    if (max(sapply(dfNames, nchar) > 8)) {
        warning("Variable names are longer than 8 characters. This could be a problem if using Mplus",
                call. = FALSE)
    }

    ## Collect basic info
    namesList <- strsplit(dfNames, split = "_")
    nSplits <- sapply(namesList, length)

    variableNames <- unique(sapply(namesList, "[[",1))
    if (length(variableNames) == 1) {
        yVar <- FALSE
    } else if (length(variableNames) == 2) {
        yVar <- TRUE
    } else if (length(variableNames) > 2) {
        stop("Variables are not named consistently or there are more than two variables. Check variable names.",
             call. = FALSE)
    }


    xNames <- namesList[sapply(namesList, "[[",1) == variableNames[[1]]]
    xWaves <- as.numeric(unique(sapply(xNames, "[[", 2)))
    xSplits <- sapply(xNames, length)
    maxWaves <- max(xWaves)
    if (yVar == TRUE) {
        yNames <- namesList[sapply(namesList, "[[", 1) == variableNames[[2]]]
        ySplits <- sapply(yNames, length)
        yWaves <- as.numeric(unique(sapply(yNames, "[[", 2)))
        maxWaves <- max(xWaves, yWaves)
    }


    ## Check whether splits make sense
    if (any(nSplits > 3)) {
        stop("Too may splits. Check variable names", call. = FALSE)
    }
    if (length(unique(sapply(xNames, length))) > 1) {
        stop("X variable names cannot be split in the same way for each variable. Check variable names.",
             call. = FALSE)
    }
    if (yVar == TRUE) {
        if (length(unique(sapply(yNames, length))) > 1) {
            stop("Y variable names cannot be split in the same way for each variable. Check variable names.",
                 call. = FALSE)
        }
    }
    if (any(nSplits == 1)) {
        stop("Not able to identify number of waves. Check variable names.", call.=FALSE)
    }
    
    

    ## Check number of indicators
    if (max(xSplits) == 3) {
        if (length(unique(table(sapply(xNames, "[[", 2)))) > 1) {
            stop("Different number of indicators per wave. Check data.", call. = FALSE)
        }
        xIndicators <- as.numeric(max(unique(sapply(xNames, "[[", 3))))
    } else {
        xIndicators <- 1
    }
    
    if (yVar == TRUE) {
        if (max(ySplits) == 3) {
            if (length(unique(table(sapply(yNames, "[[", 2)))) > 1) {
                stop("Different number of indicators per wave. Check data.", call. = FALSE)
            }
            yIndicators <- as.numeric(max(unique(sapply(yNames, "[[", 3))))
        } else {
            yIndicators <- 1
        }
    } 

    ## Create list of info
    if (yVar == TRUE) {
        info <- list(
            gen = list(
                maxWaves = maxWaves,
                yVar = TRUE
                ),
            x = list(
                name = variableNames[1],
                waves = 1:maxWaves,
                actualWaves = xWaves,
                indicators = xIndicators
            ),
            y = list(
                name = variableNames[2],
                waves = 1:maxWaves,
                actualWaves = yWaves,
                indicators = yIndicators
            )
        )
    } else {
        info <- list(
            gen = list(
                maxWaves = maxWaves,
                yVar = FALSE
                ),
            x = list(
                name = variableNames[1],
                waves = 1:maxWaves,
                actualWaves = xWaves,
                indicators = xIndicators
            ),
            y = NULL
        )
    }
    return(info)
}

#' Get summary data from lavaan output
#' @param fitObject Lavaan object
#' @noRd
.summarizeLavaan <- function(panelModel,
                            info,
                            fitObject,
                            crossLag = crossLag,
                            ma = ma,
                            clma = ma,
                            traitCors = traitCors,
                            arCors = arCors,
                            stateCors = stateCors,
                            residCors = residCors,
                            slope = slope,
                            stationarity = stationarity,
                            lags) {
    ## Calculate values for summary
    est <- lavaan::parameterEstimates(fitObject)
    est.std <- lavaan::standardizedSolution(fitObject)
    fit <- lavaan::fitMeasures(fitObject)

    ## Set Model Name
    if (panelModel == "starts") {
        mName <- "Stable Trait Autoregressive Trait State Model"
    } else if (panelModel == "riclpm") {
        mName <- "Random Intercept Cross-Lagged Panel Model"
    } else if (panelModel == "clpm") {
        mName <- "Cross-Lagged Panel Model"
    } else if (panelModel == "arts") {
        mName <- "Autoregressive Trait State Model"
    } else if (panelModel == "art") {
        mName <- "Autoregressive Trait Model"
    } else if (panelModel == "sts") {
        mName <- "Stable Trait State Model"
    } else if (panelModel == "gclm") {
        mName <- "General Cross-Lagged Model"
    } else if (panelModel == "dpm_c") {
        mName <- "Constrained Dynamic Panel Model"
    } else if (panelModel == "dpm_p") {
        mName <- "Predetermined Dynamic Panel Model"
    } else if (panelModel == "lgcm") {
        mName <- "Latent Growth Curve Model"
    } else if (panelModel == "alt") {
        mName <- "Autoregressive Latent Trajectory Model"
    } else if (panelModel == "lcmsr") {
        mName <- "Latent Curve Model with Structured Residuals"
    } else if (panelModel == "measurement") {
        mName <- "Measurement Model"
    }
    

    ## Set names
    x.name <- info$x$name
    if (info$gen$yVar == TRUE) {
        y.name <- info$y$name
    } else {
        y.name <- NULL
    }
    
    ## Variance Decomp
    if(length(est[which(est$label=="x_tVar"), "est"]) == 0) {
        trait.x <- 0
        trait.x.se <- NA
    } else {
        trait.x <- est[which(est$label=="x_tVar"), "est"]
        trait.x.se <- est[which(est$label=="x_tVar"), "se"]
    }
    if(length(est[which(est$label=="y_tVar"), "est"]) == 0) {
        trait.y <- 0
        trait.y.se <- NA
    } else {
        trait.y <- est[which(est$label=="y_tVar"), "est"]
        trait.y.se <- est[which(est$label=="y_tVar"), "se"]
    }
    if(length(est[which(est$label=="xvar1"), "est"]) == 0) {
        ar.x <- 0
        ar.x.se <- NA
    } else {
        ar.x <- est[which(est$label=="xvar1"), "est"]
        ar.x.se <- est[which(est$label=="xvar1"), "se"]
    }
    if(length(est[which(est$label=="yvar1"), "est"]) == 0) {
        ar.y <- 0
        ar.y.se <- 0
    } else {
        ar.y <- est[which(est$label=="yvar1"), "est"]
        ar.y.se <- est[which(est$label=="yvar1"), "se"]
    }
    if(length(est[which(est$label=="sx1"), "est"]) == 0) {
        state.x <- 0
        state.x.se <- NA
    } else {
        state.x <- est[which(est$label=="sx1"), "est"]
        state.x.se <- est[which(est$label=="sx1"), "se"]
    }
    if(length(est[which(est$label=="sy1"), "est"]) == 0) {
        state.y <- 0
        state.y.se <- NA
    } else {
        state.y <- est[which(est$label=="sy1"), "est"]
        state.y.se <- est[which(est$label=="sy1"), "se"]
    }
    trait.x.p <- trait.x/(sum(trait.x, ar.x, state.x))
    trait.y.p <- trait.y/(sum(trait.y, ar.y, state.y))
    ar.x.p <- ar.x/(sum(trait.x, ar.x, state.x))
    ar.y.p <- ar.y/(sum(trait.y, ar.y, state.y))
    state.x.p <- state.x/(sum(trait.x, ar.x, state.x))
    state.y.p <- state.y/(sum(trait.y, ar.y, state.y))

    var_table <- data.frame(
        Variable = c("X", "Y"),
        Trait = c(
            paste0(sprintf("%.3f", round(trait.x, 3)),
                   " (",
                   sprintf("%.3f", round(trait.x.se, 3)),
                   ")"),
            paste0(sprintf("%.3f", round(trait.y, 3)),
                   " (",
                   sprintf("%.3f", round(trait.y.se, 3)),
                   ")")
        ),
        AR = c(
            paste0(sprintf("%.3f", round(ar.x, 3)),
                   " (",
                   sprintf("%.3f", round(ar.x.se, 3)),
                   ")"),
            paste0(sprintf("%.3f", round(ar.y, 3)),
                   " (",
                   sprintf("%.3f", round(ar.y.se, 3)),
                   ")")
        ),
        State = c(
            paste0(sprintf("%.3f", round(state.x, 3)),
                   " (",
                   sprintf("%.3f", round(state.x.se, 3)),
                   ")"),
            paste0(sprintf("%.3f", round(state.y, 3)),
                   " (",
                   sprintf("%.3f", round(state.y.se, 3)),
                   ")")
        )
    )

    ## Slope
    if (length(est[which(est$label == "x_slVar"), "est"]) == 0) {
        slope.var.x <- 0
        slope.var.se.x <- NA
        slope.mean.x <- NA
        slope.mean.se.x <- NA
    } else {
        slope.var.x <- est[which(est$label == "x_slVar"), "est"]
        slope.var.se.x <- est[which(est$label == "x_slVar"), "se"]
        slope.mean.x <- est[which(est$lhs == paste0("sl_", x.name) &
            est$op == "~1"), "est"]
        slope.mean.se.x <- est[which(est$lhs == paste0("sl_", x.name) &
            est$op == "~1"), "se"]
    }

    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$label == "y_slVar"), "est"]) == 0) {
            slope.var.y <- 0
            slope.var.se.y <- NA
            slope.mean.y <- NA
            slope.mean.se.y <- NA
        } else {
            slope.var.y <- est[which(est$label == "y_slVar"), "est"]
            slope.var.se.y <- est[which(est$label == "y_slVar"), "se"]
            slope.mean.y <- est[which(est$lhs == paste0("sl_", y.name) &
                est$op == "~1"), "est"]
            slope.mean.se.y <- est[which(est$lhs == paste0("sl_", y.name) &
                est$op == "~1"), "se"]
        }
    } else {
        slope.var.y <- 0
        slope.var.se.y <- NA
        slope.mean.y <- NA
        slope.mean.se.y <- NA
    }
    


    ## Correlations
    if (info$gen$yVar == TRUE) {
        trait.cor <- est.std[which(est.std$label == "cov_txty"), "est.std"]
        trait.cor.u <- est[which(est$label == "cov_txty"), "est"]
        trait.cor.se <- est[which(est$label == "cov_txty"), "se"]
    } else {
        trait.cor <- NA 
        trait.cor.u <- NA 
        trait.cor.se <- NA 
    }

    if (info$gen$yVar == TRUE) {
        ar.cor <- est.std[which(est.std$label == "cov_ar1"), "est.std"]
        ar.cor.u <- est[which(est.std$label == "cov_ar1"), "est"]
        ar.cor.se <- est[which(est.std$label == "cov_ar1"), "se"]
    } else {
        ar.cor <- NA 
        ar.cor.u <- NA 
        ar.cor.se <- NA 
    }

    if (info$gen$yVar == TRUE) {
        state.cor <- est.std[which(est.std$label == "cov_s1"), "est.std"]
        state.cor.u <- est[which(est.std$label == "cov_s1"), "est"]
        state.cor.se <- est[which(est.std$label == "cov_s1"), "se"]
    } else {
        state.cor <- NA 
        state.cor.u <- NA 
        state.cor.se <- NA 
    }

    if (info$gen$yVar == TRUE) {
        slope.cor <- est.std[which(est.std$label == "cov_sxsy"), "est.std"]
        slope.cor.u <- est[which(est.std$label == "cov_sxsy"), "est"]
        slope.cor.se <- est[which(est.std$label == "cov_sxsy"), "se"]
    } else {
        slope.cor <- NA 
        slope.cor.u <- NA 
        slope.cor.se <- NA 
    }
    
    ## Stability
    if (length(est[which(est$label == "a_1_2"), "est"]) > 0) {
        x.stab <- est.std[which(est.std$label == "a_1_2"), "est.std"]
        x.stab.u <- est[which(est$label == "a_1_2"), "est"]
        x.stab.se <- est[which(est$label == "a_1_2"), "se"]
    } else {
        x.stab <- NA
        x.stab.u <- NA
        x.stab.se <- NA
    }
    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$label == "b_1_2"), "est"]) > 0) {
            y.stab <- est.std[which(est.std$label == "b_1_2"), "est.std"]
            y.stab.u <- est[which(est$label == "b_1_2"), "est"]
            y.stab.se <- est[which(est$label == "b_1_2"), "se"]
        } else {
            y.stab <- NA
            y.stab.u <- NA
            y.stab.se <- NA
        }
    } else {
        y.stab <- NA
        y.stab.u <- NA
        y.stab.se <- NA
    }
    
    ## Cross-lags
    if (length(est[which(est$label == "c_1_2"), "est"]) > 0) {
        yOnX <- est.std[which(est.std$label == "c_1_2"), "est.std"]
        yOnX.u <- est[which(est$label == "c_1_2"), "est"]
        yOnX.se <- est[which(est$label == "c_1_2"), "se"]
    } else {
        yOnX <- NA
        yOnX.u <- NA
        yOnX.se <- NA
    }
    if (length(est[which(est$label == "d_1_2"), "est"]) > 0) {
        xOnY <- est.std[which(est.std$label == "d_1_2"), "est.std"]
        xOnY.u <- est[which(est$label == "d_1_2"), "est"]
        xOnY.se <- est[which(est$label == "d_1_2"), "se"]
    } else {
        xOnY <- NA
        xOnY.u <- NA
        xOnY.se <- NA
    }

    ## Moving Averages
    if (length(est[which(est$label == "ma_a2"), "est"]) > 0) {
        x_ma <- est.std[which(est.std$label == "ma_a2"), "est.std"]
        x_ma.u <- est[which(est$label == "ma_a2"), "est"]
        x_ma.se <- est[which(est$label == "ma_a2"), "se"]
    } else {
        x_ma <- NA
        x_ma.u <- NA
        x_ma.se <- NA
    }
    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$label == "ma_b2"), "est"]) > 0) {
            y_ma <- est.std[which(est.std$label == "ma_b2"), "est.std"]
            y_ma.u <- est[which(est$label == "ma_b2"), "est"]
            y_ma.se <- est[which(est$label == "ma_b2"), "se"]
        } else {
            y_ma <- NA
            y_ma.u <- NA
            y_ma.se <- NA
        }
    } else {
        y_ma <- NA
        y_ma.u <- NA
        y_ma.se <- NA
    }

    ## Cross-lag Moving Averages
    if (length(est[which(est$label == "ma_c2"), "est"]) > 0) {
        yOnX_ma <- est.std[which(est.std$label == "ma_c2"), "est.std"]
        yOnX_ma.u <- est[which(est$label == "ma_c2"), "est"]
        yOnX_ma.se <- est[which(est$label == "ma_c2"), "se"]
    } else {
        yOnX_ma <- NA
        yOnX_ma.u <- NA
        yOnX_ma.se <- NA
    }
    if (length(est[which(est$label == "ma_d2"), "est"]) > 0) {
        xOnY_ma <- est.std[which(est.std$label == "ma_d2"), "est.std"]
        xOnY_ma.u <- est[which(est$label == "ma_d2"), "est"]
        xOnY_ma.se <- est[which(est$label == "ma_d2"), "se"]
    } else {
        xOnY_ma <- NA
        xOnY_ma.u <- NA
        xOnY_ma.se <- NA
    }
    
    ## Fit
    chi2 <- fit[["chisq"]]
    chi2df <- fit[["df"]]
    chi2p <- fit[["pvalue"]]
    cfi <- fit[["cfi"]]
    tli <- fit[["tli"]]
    srmr <- fit[["srmr"]]
    rmsea <- fit[["rmsea"]]
    aic <- fit[["aic"]]
    bic <- fit[["bic"]]

    
    outputList <- list(
        info = info,
        model = panelModel,
        mName = mName,
        program = "Lavaan",
        crossLag = crossLag,
        lags = lags,
        ma = ma,
        clma = clma,
        traitCors = traitCors,
        arCors = arCors,
        stateCors = stateCors,
        residCors = residCors,
        stationarity = stationarity,
        slope = slope,
        x.name = x.name,
        y.name = y.name,
        trait.x = trait.x.p,
        trait.y = trait.y.p,
        ar.x = ar.x.p,
        ar.y = ar.y.p,
        state.x = state.x.p,
        state.y = state.y.p,
        var_table = var_table,
        slope.var.x = slope.var.x,
        slope.var.se.x = slope.var.se.x,
        slope.mean.x = slope.mean.x,
        slope.mean.se.x = slope.mean.se.x,
        slope.var.y = slope.var.y,
        slope.var.se.y = slope.var.se.y,
        slope.mean.y = slope.mean.y,
        slope.mean.se.y = slope.mean.se.y,
        trait.cor = trait.cor,
        trait.cor.u = trait.cor.u,
        trait.cor.se = trait.cor.se,
        slope.cor = slope.cor,
        slope.cor.u = slope.cor.u,
        slope.cor.se = slope.cor.se,
        ar.cor = ar.cor,
        ar.cor.u = ar.cor.u,
        ar.cor.se = ar.cor.se,
        state.cor = state.cor,
        state.cor.u = state.cor.u,
        state.cor.se = state.cor.se,
        x.stab = x.stab,
        x.stab.u = x.stab.u,
        x.stab.se = x.stab.se,
        y.stab = y.stab,
        y.stab.u = y.stab.u,
        y.stab.se = y.stab.se,
        yOnX = yOnX,
        xOnY = xOnY,
        yOnX.u = yOnX.u,
        xOnY.u = xOnY.u,
        yOnX.se = yOnX.se,
        xOnY.se = xOnY.se,
        x_ma = x_ma,
        x_ma.u = x_ma.u,
        x_ma.se = x_ma.se,
        y_ma = y_ma,
        y_ma.u = y_ma.u,
        y_ma.se = y_ma.se,
        yOnX_ma = yOnX_ma,
        yOnX_ma.u = yOnX_ma.u,
        yOnX_ma.se = yOnX_ma.se,
        xOnY_ma = xOnY_ma,
        xOnY_ma.u = xOnY_ma.u,
        xOnY_ma.se = xOnY_ma.se,
        chi2 = chi2,
        chi2df = chi2df,
        chi2p = chi2p,
        cfi = cfi,
        tli = tli,
        srmr = srmr,
        rmsea = rmsea,
        aic = aic,
        bic = bic)
    class(outputList) <- "pcSum"
    return(outputList)
}

#' Summarize info from mplus object
#' @param info information object generated from `getInfo()`
#' @param fitObject mplus object with output
#' @noRd
.summarizeMplus <- function(panelModel,
                            info,
                            fitObject,
                            crossLag = crossLag,
                            ma = ma,
                            clma = ma,
                            traitCors = traitCors,
                            arCors = arCors,
                            stateCors = stateCors,
                            residCors = residCors,
                            slope = slope,
                            stationarity = stationarity,
                            lags
                            ) {
    ## Extract estimates and fit info
    est <- fitObject$results$parameters$unstandardized
    est.std <- fitObject$results$parameters$stdyx.standardized
    fit <- fitObject$results$summaries

    ## Set Model Name
    if (panelModel == "starts") {
        mName <- "Stable Trait Autoregressive Trait State Model"
    } else if (panelModel == "riclpm") {
        mName <- "Random Intercept Cross-Lagged Panel Model"
    } else if (panelModel == "clpm") {
        mName <- "Cross-Lagged Panel Model"
    } else if (panelModel == "arts") {
        mName <- "Autoregressive Trait State Model"
    } else if (panelModel == "art") {
        mName <- "Autoregressive Trait Model"
    } else if (panelModel == "sts") {
        mName <- "Stable Trait State Model"
    } else if (panelModel == "gclm") {
        mName <- "General Cross-Lagged Model"
    } else if (panelModel == "dpm_c") {
        mName <- "Constrained Dynamic Panel Model"
    } else if (panelModel == "dpm_p") {
        mName <- "Predetermined Dynamic Panel Model"
    } else if (panelModel == "lgcm") {
        mName <- "Latent Growth Curve lModel"
    } else if (panelModel == "alt") {
        mName <- "Autoregressive Latent Trajectory Model"
    } else if (panelModel == "lcmsr") {
        mName <- "Latent Curve Model with Structured Residuals Model"
    } else if (panelModel == "measurement") {
        mName <- "Measurement Model"
    }
    
    



    ## Set names
    x.name <- info$x$name
    if (info$gen$yVar == TRUE) {
        y.name <- info$y$name
    } else {
        y.name <- NULL
    }
    
    t.x  <- paste("T", toupper(info$x$name), sep = "_")
    if (info$gen$yVar == TRUE) {
        t.y <- paste("T", toupper(info$y$name), sep = "_")
    } 
    ## Slope
    sl.x <- paste("SL", toupper(info$x$name), sep = "_")
    if (info$gen$yVar == TRUE) {
        sl.y <- paste("SL", toupper(info$y$name), sep = "_")
    }
    ## AR Variance is captured in impulses
    a.x  <- paste("A", toupper(info$x$name), "1", sep = "_")
    a.x2  <- paste("A", toupper(info$x$name), "2", sep = "_")
    if (info$gen$yVar == TRUE) {
        a.y <- paste("A", toupper(info$y$name), "1", sep = "_")
        a.y2 <- paste("A", toupper(info$y$name), "2", sep = "_")
    }
    i.x  <- paste("I", toupper(info$x$name), "1", sep = "_")
    i.x2  <- paste("I", toupper(info$x$name), "2", sep = "_")
    if (info$gen$yVar == TRUE) {
        i.y <- paste("I", toupper(info$y$name), "1", sep = "_")
        i.y2 <- paste("I", toupper(info$y$name), "2", sep = "_")
    }
    s.x  <- paste("S", toupper(info$x$name), "1", sep = "_")
    if (info$gen$yVar == TRUE) {
        s.y <- paste("S", toupper(info$y$name), "1", sep = "_")
    }
    
    ## Variance Decomp
    if (length(est[which(est$paramHeader == "Variances" &
                         est$param == t.x), "est"]) == 0) {
        trait.x <- 0
        trait.x.se <- NA
    } else {
        trait.x <- est[which(est$paramHeader == "Variances" &
                             est$param == t.x), "est"]
        trait.x.se <- est[which(est$paramHeader == "Variances" &
                             est$param == t.x), "se"]
    }
    
    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$paramHeader == "Variances" &
                             est$param == t.y), "est"]) == 0) {
            trait.y <- 0
            trait.y.se <- NA
        } else {
            trait.y <- est[which(est$paramHeader == "Variances" &
                                 est$param == t.y), "est"]
            trait.y.se <- est[which(est$paramHeader == "Variances" &
                                 est$param == t.y), "se"]
        }
    } else {
        trait.y <- NA
        trait.y.se <- NA
    }

    if (length(est[which(est$paramHeader == "Variances" &
                         est$param == i.x), "est"]) == 0) {
        ar.x <- 0
        ar.x.se <- NA
    } else {
        ar.x <- est[which(est$paramHeader == "Variances" &
                          est$param == i.x), "est"]
        ar.x.se <- est[which(est$paramHeader == "Variances" &
                          est$param == i.x), "se"]
    }
    
    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$paramHeader == "Variances" &
                             est$param == i.y), "est"]) == 0) {
            ar.y <- 0
            ar.y.se <- NA
        } else {
            ar.y <- est[which(est$paramHeader == "Variances" &
                              est$param == i.y), "est"]
            ar.y.se <- est[which(est$paramHeader == "Variances" &
                              est$param == i.y), "se"]
        }
    } else {
        ar.y <- NA
        ar.y.se <- NA
    }
    if (length(est[which(est$paramHeader == "Variances" &
                         est$param == s.x), "est"]) == 0) {
        state.x <- 0
        state.x.se <- NA
    } else {
        state.x <- est[which(est$paramHeader == "Variances" &
                             est$param == s.x), "est"]
        state.x.se <- est[which(est$paramHeader == "Variances" &
                             est$param == s.x), "se"]
    }
    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$paramHeader == "Variances" &
                             est$param == s.y), "est"]) == 0) {
            state.y <- 0
            state.y.se <- 0
        } else {
            state.y <- est[which(est$paramHeader == "Variances" &
                                 est$param == s.y), "est"]
            state.y.se <- est[which(est$paramHeader == "Variances" &
                                 est$param == s.y), "se"]
        }
    } else {
        state.y <- NA
        state.y.se <- NA
    }
    trait.x.p <- trait.x/(sum(trait.x, ar.x, state.x))
    trait.y.p <- trait.y/(sum(trait.y, ar.y, state.y))
    ar.x.p <- ar.x/(sum(trait.x, ar.x, state.x))
    ar.y.p <- ar.y/(sum(trait.y, ar.y, state.y))
    state.x.p <- state.x/(sum(trait.x, ar.x, state.x))
    state.y.p <- state.y/(sum(trait.y, ar.y, state.y))

    var_table <- data.frame(
        Variable = c("X", "Y"),
        Trait = c(
            paste0(sprintf("%.3f", round(trait.x, 3)),
                   " (",
                   sprintf("%.3f", round(trait.x.se, 3)),
                   ")"),
            paste0(sprintf("%.3f", round(trait.y, 3)),
                   " (",
                   sprintf("%.3f", round(trait.y.se, 3)),
                   ")")
        ),
        AR = c(
            paste0(sprintf("%.3f", round(ar.x, 3)),
                   " (",
                   sprintf("%.3f", round(ar.x.se, 3)),
                   ")"),
            paste0(sprintf("%.3f", round(ar.y, 3)),
                   " (",
                   sprintf("%.3f", round(ar.y.se, 3)),
                   ")")
        ),
        State = c(
            paste0(sprintf("%.3f", round(state.x, 3)),
                   " (",
                   sprintf("%.3f", round(state.x.se, 3)),
                   ")"),
            paste0(sprintf("%.3f", round(state.y, 3)),
                   " (",
                   sprintf("%.3f", round(state.y.se, 3)),
                   ")")
        )
    )
    

    ## Slope
    if (length(est[which(est$paramHeader == "Variances" &
                         est$param == sl.x), "est"]) == 0) {
        slope.var.x <- 0
        slope.var.se.x <- NA
    } else {
        slope.var.x <- est[which(est$paramHeader == "Variances" &
                                 est$param == sl.x), "est"]
        slope.var.se.x <- est[which(est$paramHeader == "Variances" &
                             est$param == sl.x), "se"]
    }
    if (length(est[which(est$paramHeader == "Means" &
                         est$param == sl.x), "est"]) == 0) {
        slope.mean.x <- 0
        slope.mean.se.x <- NA
    } else {
        slope.mean.x <- est[which(est$paramHeader == "Means" &
                                  est$param == sl.x), "est"]
        slope.mean.se.x <- est[which(est$paramHeader == "Means" &
                             est$param == sl.x), "se"]
    }
    if (info$gen$yVar == TRUE) {
        if (length(est[which(est$paramHeader == "Variances" &
            est$param == sl.y), "est"]) == 0) {
            slope.var.y <- 0
            slope.var.se.y <- NA
        } else {
            slope.var.y <- est[which(est$paramHeader == "Variances" &
                                     est$param == sl.y), "est"]
            slope.var.se.y <- est[which(est$paramHeader == "Variances" &
                                     est$param == sl.y), "se"]
        }
        if (length(est[which(est$paramHeader == "Means" &
            est$param == sl.y), "est"]) == 0) {
            slope.mean.y <- 0
            slope.mean.se.y <- NA
        } else {
            slope.mean.y <- est[which(est$paramHeader == "Means" &
                                      est$param == sl.y), "est"]
            slope.mean.se.y <- est[which(est$paramHeader == "Means" &
                                      est$param == sl.y), "se"]
        }
    } else {
        slope.var.y <- NULL
        slope.var.se.y <- NULL
        slope.mean.y <- NULL
        slope.mean.se.y <- NULL
    }
         
    
    ## Correlations
    if (info$gen$yVar == TRUE) {
        trait.cor <- est.std[which(est.std$paramHeader == paste0(t.x, ".WITH") &
                                   est.std$param == t.y), "est"]
        trait.cor.u <- est[which(est$paramHeader == paste0(t.x, ".WITH") &
                                   est$param == t.y), "est"]
        trait.cor.se <- est[which(est$paramHeader == paste0(t.x, ".WITH") &
                                   est$param == t.y), "se"]
    } else {
        trait.cor <- NA
        trait.cor.u <- NA
        trait.cor.se <- NA
    }

    if (info$gen$yVar == TRUE) {
        slope.cor <- est.std[which(est.std$paramHeader == paste0(sl.x, ".WITH") &
                                   est.std$param == sl.y), "est"]
        slope.cor.u <- est[which(est$paramHeader == paste0(sl.x, ".WITH") &
                                   est$param == sl.y), "est"]
        slope.cor.se <- est[which(est$paramHeader == paste0(sl.x, ".WITH") &
                                   est$param == sl.y), "se"]
    } else {
        slope.cor <- NA
        slope.cor.u <- NA
        slope.cor.se <- NA
    }
    
    if (info$gen$yVar == TRUE) {
        ar.cor <- est.std[which(est.std$paramHeader == paste0(i.x, ".WITH") &
                                est.std$param == i.y), "est"]
        ar.cor.u <- est[which(est$paramHeader == paste0(i.x, ".WITH") &
                                est$param == i.y), "est"]
        ar.cor.se <- est[which(est$paramHeader == paste0(i.x, ".WITH") &
                                est$param == i.y), "se"]
    } else {
        ar.cor <- NA
        ar.cor.u <- NA
        ar.cor.se <- NA
    }
    if (info$gen$yVar == TRUE) {
        state.cor <- est.std[which(est.std$paramHeader == paste0(s.x, ".WITH") &
                                   est.std$param == s.y), "est"]
        state.cor.u <- est[which(est$paramHeader == paste0(s.x, ".WITH") &
                                   est$param == s.y), "est"]
        state.cor.se <- est[which(est$paramHeader == paste0(s.x, ".WITH") &
                                   est$param == s.y), "se"]
    } else {
        state.cor <- NA
        state.cor.u <- NA
        state.cor.se <- NA
    }
    
    ## Stability
    if (length(est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                             est.std$param == a.x), "est"]) == 0) {
        x.stab <- NA
        x.stab.u <- NA
        x.stab.se <- NA
    } else {
        x.stab <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                                est.std$param == a.x), "est"]
        x.stab.u <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                                est$param == a.x), "est"]
        x.stab.se <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                                est$param == a.x), "se"]
    }
    if (info$gen$yVar == TRUE) {
        if (length(est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                 est.std$param == a.y), "est"]) == 0) {
            y.stab <- NA
            y.stab.u <- NA
            y.stab.se <- NA
        } else {
            y.stab <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                    est.std$param == a.y), "est"]
            y.stab.u <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                    est$param == a.y), "est"]
            y.stab.se <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                    est$param == a.y), "se"]
        }
    } else {
        y.stab <- NA
        y.stab.u <- NA
        y.stab.se <- NA
    }
    
    ## Cross-lags
    if (info$gen$yVar == TRUE) {
        if (length(est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                 est.std$param == a.x), "est"]) > 0 &
            length(est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                                 est.std$param == a.y), "est"]) > 0) {
            yOnX <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                  est.std$param == a.x), "est"]
            xOnY <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                                  est.std$param == a.y), "est"]
            yOnX.u <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                  est$param == a.x), "est"]
            xOnY.u <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                                  est$param == a.y), "est"]
            yOnX.se <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                  est$param == a.x), "se"]
            xOnY.se <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                                  est$param == a.y), "se"]
        } else {
            yOnX <- NA
            xOnY <- NA
            yOnX.u <- NA
            xOnY.u <- NA
            yOnX.se <- NA
            xOnY.se <- NA
        }
    } else {
        yOnX <- NA
        xOnY <- NA
        yOnX.u <- NA
        xOnY.u <- NA
        yOnX.se <- NA
        xOnY.se <- NA
    }

    ## Moving Averages
    if (length(est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                             est.std$param == i.x), "est"]) > 0) {
        x_ma <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                              est.std$param == i.x), "est"]
        x_ma.u <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                            est$param == i.x), "est"]
        x_ma.se <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                             est$param == i.x), "se"]
    } else {
        x_ma <- NA
        x_ma.u <- NA
        x_ma.se <- NA
    }

    if (info$gen$yVar == TRUE) {
        if (length(est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                             est.std$param == i.y), "est"]) > 0) {
        y_ma <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                 est.std$param == i.y), "est"]
        y_ma.u <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                               est$param == i.y), "est"]
        y_ma.se <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                est$param == i.y), "se"]
        } else {
            y_ma <- NA
            y_ma.u <- NA
            y_ma.se <- NA
        }
    } else {
        y_ma <- NA
        y_ma.u <- NA
        y_ma.se <- NA
    }
    
    


    
    ## Cross-Lagged Moving Averages
    if (info$gen$yVar == TRUE) {
        if (length(est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                 est.std$param == i.x), "est"]) > 0 &
            length(est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                                 est.std$param == i.y), "est"]) > 0) {
            yOnX_ma <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                  est.std$param == i.x), "est"]
            xOnY_ma <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                                  est.std$param == i.y), "est"]
            yOnX_ma.u <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                  est$param == i.x), "est"]
            xOnY_ma.u <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                                  est$param == i.y), "est"]
            yOnX_ma.se <- est[which(est$paramHeader == paste0(a.y2, ".ON") &
                                  est$param == i.x), "se"]
            xOnY_ma.se <- est[which(est$paramHeader == paste0(a.x2, ".ON") &
                                  est$param == i.y), "se"]
        } else {
            yOnX_ma <- NA
            xOnY_ma <- NA
            yOnX_ma.u <- NA
            xOnY_ma.u <- NA
            yOnX_ma.se <- NA
            xOnY_ma.se <- NA
        }
    } else {
        yOnX_ma <- NA
        xOnY_ma <- NA
        yOnX_ma.u <- NA
        xOnY_ma.u <- NA
        yOnX_ma.se <- NA
        xOnY_ma.se <- NA
    }
    
    ## Fit
    chi2 <- fit[["ChiSqM_Value"]]
    chi2df <- fit[["ChiSqM_DF"]]
    chi2p <- fit[["ChiSqM_PValue"]]
    cfi <- fit[["CFI"]]
    tli <- fit[["TLI"]]
    srmr <- fit[["SRMR"]]
    rmsea <- fit[["RMSEA_Estimate"]]
    aic <- fit[["AIC"]]
    bic <- fit[["BIC"]]
    outputList <- list(
        info = info,
        model = panelModel,
        mName = mName,
        program = "Mplus",
        crossLag = crossLag,
        lags = lags,
        ma = ma,
        clma = clma,
        traitCors = traitCors,
        arCors = arCors,
        stateCors = stateCors,
        residCors = residCors,
        stationarity = stationarity,
        slope = slope,
        x.name = x.name,
        y.name = y.name,
        trait.x = trait.x.p,
        trait.y = trait.y.p,
        ar.x = ar.x.p,
        ar.y = ar.y.p,
        state.x = state.x.p,
        state.y = state.y.p,
        var_table = var_table,
        slope.var.x = slope.var.x,
        slope.var.se.x = slope.var.se.x,
        slope.mean.x = slope.mean.x,
        slope.mean.se.x = slope.mean.se.x,
        slope.var.y = slope.var.y,
        slope.var.se.y = slope.var.se.y,
        slope.mean.y = slope.mean.y,
        slope.mean.se.y = slope.mean.se.y,
        trait.cor = trait.cor,
        trait.cor.u = trait.cor.u,
        trait.cor.se = trait.cor.se,
        slope.cor = slope.cor,
        slope.cor.u = slope.cor.u,
        slope.cor.se = slope.cor.se,
        ar.cor = ar.cor,
        ar.cor.u = ar.cor.u,
        ar.cor.se = ar.cor.se,
        state.cor = state.cor,
        state.cor.u = state.cor.u,
        state.cor.se = state.cor.se,
        x.stab = x.stab,
        x.stab.u = x.stab.u,
        x.stab.se = x.stab.se,
        y.stab = y.stab,
        y.stab.u = y.stab.u,
        y.stab.se = y.stab.se,
        yOnX = yOnX,
        xOnY = xOnY,
        yOnX.u = yOnX.u,
        xOnY.u = xOnY.u,
        yOnX.se = yOnX.se,
        xOnY.se = xOnY.se,
        x_ma = x_ma,
        x_ma.u = x_ma.u,
        x_ma.se = x_ma.se,
        y_ma = y_ma,
        y_ma.u = y_ma.u,
        y_ma.se = y_ma.se,
        yOnX_ma = yOnX_ma,
        yOnX_ma.u = yOnX_ma.u,
        yOnX_ma.se = yOnX_ma.se,
        xOnY_ma = xOnY_ma,
        xOnY_ma.u = xOnY_ma.u,
        xOnY_ma.se = xOnY_ma.se,
        chi2 = chi2,
        chi2df = chi2df,
        chi2p = chi2p,
        cfi = cfi,
        tli = tli,
        srmr = srmr,
        rmsea = rmsea,
        aic = aic,
        bic = bic)
    class(outputList) <- "pcSum"
    return(outputList)
}


#' Prints results from `panelcoder()`
#'
#' @param x Results from `panelcoder()`.
#' @param ... Additional arguments to `print()`.
#'
#' @export
print.pcSum <- function(x, ...) {
    cat(x$mName, "\n")
    cat("Estimated with ", x$program, "\n")
    cat("X Variable: ", x$x.name, "\n")
    if (x$info$gen$yVar == TRUE) {
        cat("Y Variable: ", x$y.name, "\n")
    }
    cat("Number of Waves: ", x$info$gen$maxWaves, "\n")
    cat("\n")
    cat("Model Summary: \n")
    cat("Model: Chi2 (df = ",
        sprintf("%i", x$chi2df), ") = ",
        sprintf("%.3f", x$chi2), ", p = ",
        sprintf("%.3f", x$chi2p), "\n",
        sep="")
    cat("\n")
    cat("Fit Indices:")
    cat("\n")
    cat("CFI = ",
        sprintf("%.3f", x$cfi),
        ", TLI = ",
        sprintf("%.3f", x$tli),
        ", SRMR = ",
        sprintf("%.3f", x$srmr), "\n",
        sep="")
    cat("RMSEA = ",
        sprintf("%.3f", x$rmsea),
        "\n")
    cat("AIC = ",
        sprintf("%.3f", x$aic),
        ", BIC = ",
        sprintf("%.3f", x$bic), "\n")
    cat("\n")
    cat("\n")
    cat("Estimates (First Wave if No Stationarity) \n")
    cat("\n")
    cat("Variances: \n")
    if (x$info$gen$yVar == TRUE){
        print(x$var_table)
    } else {
        print(x$var_table[1,])
    }
    cat("\n")
    if (x$stationarity == "full") {
        cat("Variance Decompostion (With Stationarity) \n")
        cat("X Variable: \n")
        cat("Trait: ",
            sprintf("%.2f", 100*x$trait.x),
            "%, AR: ",
            sprintf("%.2f", 100*x$ar.x),
            "%, State: ",
            sprintf("%.2f", 100*x$state.x),
            "% \n",
            sep="")
        if (x$info$gen$yVar == TRUE) {
            cat("Y Variable: \n")
            cat("Trait: ",
                sprintf("%.2f", 100*x$trait.y),
                "%, AR: ",
                sprintf("%.2f", 100*x$ar.y),
                "%, State: ",
                sprintf("%.2f", 100*x$state.y),
                "% \n",
                sep="")
        }
        cat("\n")
    }
    if (x$slope != "none") {
        cat("Slopes: \n")
        cat("X:   ")
        cat("Variance: ", sprintf("%.3f", x$slope.var.x), " (",
            sprintf("%.3f", x$slope.var.se.x), 
            ")  Mean: ",
            sprintf("%.3f", x$slope.mean.x), " (",
            sprintf("%.3f", x$slope.mean.se.x),
            ") \n", sep = "")
        if (x$info$gen$yVar == TRUE) {
            cat("Y:   ")
            cat("Variance: ", sprintf("%.3f", x$slope.var.y), " (",
                sprintf("%.3f", x$slope.var.se.y), 
                ")  Mean: ",
                sprintf("%.3f", x$slope.mean.y), " (",
                sprintf("%.3f", x$slope.mean.se.y),
                ") \n", sep = "")
        }
        cat("\n")
    }
    if (x$model != "lgcm") {
        cat("Stability: \n")
        cat("X: ", sprintf("%.3f", x$x.stab.u), " (",
            sprintf("%.3f", x$x.stab.se), "), ",
            "Standardized: ", sprintf("%.3f", x$x.stab),
            "\n", sep = "")
        if (x$info$gen$yVar == TRUE) {
            cat("Y: ", sprintf("%.3f", x$y.stab.u), " (",
                sprintf("%.3f", x$y.stab.se), "), ",
                "Standardized: ", sprintf("%.3f", x$y.stab),
                "\n", sep = "")
        }
        cat("\n")
    }
    if (!is.na(x$yOnX) & !is.na(x$xOnY)) {
        cat("Cross-Lag Paths: \n")
        cat("Y predicted from X: ",
            sprintf("%.3f", round(x$yOnX.u, 3)),
            " (", sprintf("%.3f", round(x$yOnX.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$yOnX, 3)),
            "\n",
            sep = ""
        )
        cat("X predicted from Y: ",
            sprintf("%.3f", round(x$xOnY.u, 3)), " (",
            sprintf("%.3f", round(x$xOnY.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$xOnY, 3)),
            "\n",
            sep = ""
            )
        cat("\n")
    }

    if (!is.na(x$x_ma)) {
        cat("Moving Averages: \n")
        cat("X: ",
            sprintf("%.3f", round(x$x_ma.u, 3)), " (",
            sprintf("%.3f", round(x$x_ma.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$x_ma, 3)),
            "\n",
            sep = ""
            )
        if (x$info$gen$yVar == TRUE) {
            cat("Y: ",
            sprintf("%.3f", round(x$y_ma.u, 3)), " (",
            sprintf("%.3f", round(x$y_ma.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$y_ma, 3)),
            "\n",
            sep = ""
            )
        }
        cat("\n")
    }

    if (!is.na(x$xOnY_ma)) {
        cat("Cross-Lagged Moving Averages: \n")
        cat("X Predicted from Y: ",
            sprintf("%.3f", round(x$xOnY_ma.u, 3)), " (",
            sprintf("%.3f", round(x$xOnY_ma.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$xOnY_ma, 3)),
            "\n",
            sep = ""
            )
        cat("Y Predicted from X: ",
            sprintf("%.3f", round(x$yOnX_ma.u, 3)), " (",
            sprintf("%.3f", round(x$yOnX_ma.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$yOnX_ma, 3)),
            "\n",
            sep = ""
            )
        cat("\n")
    }

    
    cat("Correlations: \n")
    cat(
        "Stable Trait: ",
        sprintf("%.3f", round(x$trait.cor.u, 3)), " (",
        sprintf("%.3f", round(x$trait.cor.se, 3)), "), Standardized: ",
        sprintf("%.3f", round(x$trait.cor, 3)),
        "\n",
        sep = ""
    )
    cat("Autoregressive Trait: ",
        sprintf("%.3f", round(x$ar.cor.u, 3)), " (",
        sprintf("%.3f", round(x$ar.cor.se, 3)), "), Standardized: ",
        sprintf("%.3f", round(x$ar.cor, 3)),
        "\n",
        sep = ""
    )
    cat(
        "State: ",
        sprintf("%.3f", round(x$state.cor.u, 3)), " (",
        sprintf("%.3f", round(x$state.cor.se, 3)), "), Standardized: ",
        sprintf("%.3f", round(x$state.cor, 3)),
        "\n",
        sep = ""
    )
    if (x$slope != "none") {
        cat(
            "Slope: ",
            sprintf("%.3f", round(x$slope.cor.u, 3)), " (",
            sprintf("%.3f", round(x$slope.cor.se, 3)), "), Standardized: ",
            sprintf("%.3f", round(x$slope.cor, 3)),
            "\n",
            sep = ""
        )
    }
    cat("\n")
    cat("\n")
}





getObservedCors <- function(df, info, var) {
    waves <- info$gen$maxWaves
    aWaves <- info[[var]]$actualWaves
    phantom <- setdiff(info[[var]]$waves, info[[var]]$actualWaves)
    name <- info[[var]]$name
    oldData <- df[, paste(name, aWaves, sep = "_")]
    if (length(phantom) > 0) {
        newData <- as.data.frame(matrix(NA,
                                        nrow = nrow(df),
                                        ncol = length(phantom)))
        names(newData) <- paste(name, phantom, sep = "_")
        finalData <- cbind(oldData, newData)
    } else {
        finalData <- oldData
    }
    corMat <- cor(finalData[,paste(name, 1:waves, sep="_")],
                  use = "pair")
    return(corMat)
}




getLavaanImpliedCors <- function(lavFit, vars, waves) {
    xNames <- paste(rep("l", waves),
                    vars[[1]],
                    1:waves,
                    sep = "_")
    if (length(vars) == 2) {
        yNames <- paste(rep("l", waves),
                        vars[[2]],
                        1:waves,
                        sep = "_")
    }
    allCors <- lavaan::lavInspect(lavFit, what = "cor.lv")
    xCors <- allCors[xNames, xNames]
    xSum <- summarizeR(xCors)
    if (length(vars) == 2) {
        yCors <- allCors[yNames, yNames]
        ySum <- summarizeR(yCors)
    }
    if (length(vars) == 1) {
        return(xSum)
    } else if (length(vars == 2)) {
        return(cbind(xSum, ySum))
    } else {
        stop("No more than two variables allowed")
    }
    return(lCors)
}

getMplusImpliedCors <- function(latCorEst, vars, waves) {
    xNames <- paste(rep("L", waves),
                    vars[[1]],
                    1:waves,
                    sep = "_")
    if (length(vars) == 2) {
        yNames <- paste(rep("L", waves),
                        vars[[2]],
                        1:waves,
                        sep = "_")
    } 
    latCorEst <- as.matrix(Matrix::forceSymmetric(latCorEst, uplo = "L"))
    xCors <- latCorEst[xNames, xNames]
    xSum <- summarizeR(xCors)
    if (length(vars) == 2) {
        yCors <- latCorEst[yNames, yNames]
        ySum <- summarizeR(yCors)
    }
    if (length(vars) == 1) {
        return(xSum)
    } else if (length(vars == 2)) {
        return(cbind(xSum, ySum))
    } else {
        stop("No more than two variables allowed")
    }
}




summarizeR <- function(corMat, nvars=1) {

    averageRs <- matrix(nrow=(nrow(corMat)/nvars-1),
                        ncol=nvars)

    for (k in 1:nvars) {
        for (i in 1:((nrow(corMat)/nvars)-1)) {
            sumR <- 0
            nValid <- 0
            for (j in seq(1, (nrow(corMat)-nvars*i), by=nvars)) {
                if(!is.na(corMat[(j+(i*nvars)+(k-1)), j+(k-1)])) {
                    sumR <- sumR + corMat[(j+(i*nvars)+(k-1)), j+(k-1)]
                    nValid <- nValid + 1
                }
                
            }
            averageRs[i,k] <- sumR/(nValid)
        }
    }
    return(averageRs)
}


combineCors <- function(df, info, program, fitObject, latent) {
    if (info$gen$yVar == TRUE) {
        varNames <- c(info$x$name, info$y$name)
        vars <- c("x", "y")
    } else {
        varNames <- c(info$x$name)
        vars <- "x"
    }
    
    ## Get implied correlations
    if (program == "mplus") {
        impliedCors <- getMplusImpliedCors(fitObject$results$tech4$latCorEst,
                                       toupper(varNames),
                                       info$gen$maxWaves)
    } else {
        impliedCors <- getLavaanImpliedCors(fitObject,
                                            varNames,
                                            info$gen$maxWaves)
    }

    impliedCors <- as.data.frame(impliedCors)
    names(impliedCors) = paste("i", vars, sep="_")

    allCors <- impliedCors

    ## Get observed correlations

    if (latent==FALSE) {
        if (length(vars) == 1) {
            observedCors <- summarizeR(getObservedCors(df, info, "x"))
        } else {
            observedCors <- cbind(summarizeR(getObservedCors(df, info, "x")),
                                  summarizeR(getObservedCors(df, info, "y")))
        }
        observedCors <- as.data.frame(observedCors)
        names(observedCors) = paste("o", vars, sep="_")
        allCors <- cbind(impliedCors, observedCors)
    }
    return(allCors)
}

    

plotCors <- function(cors) {
    lags <- nrow(cors)
    cors <- cors %>%
        dplyr::mutate(lag=dplyr::row_number()) %>%
        tidyr::pivot_longer(!lag) %>%
        tidyr::separate_wider_delim(name, " ", names = c("Source", "Variable"))
    minCor <- min(cors$value)

    ggplot2::ggplot(ggplot2::aes(x=lag,
                                 y=value,
                                 linetype=Source,
                                 color=Variable),
                    data=cors) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::ylim(min(minCor,0), 1)
}

#' Plots implied and actual correlations
#'
#' `panelPlot()` plots implied and actual stability coefficients for
#' increasingly long lags for the variables analyzed by the panelcoder
#' function. The function takes a panelcoder output object as an argument and
#' then plots these stabilities. If multiple indicators are used in the model,
#' then it will necessary to pass the results of a measurement model to compare
#' observed stabilities with implied stabilities.
#'
#' @importFrom methods is
#' @param pcOutput An object created from running the `panelcoder()` command.
#' @param measurement An object created from running the `panelcoder()` command
#'   specifying a measurement model.
#' @export
panelPlot <- function(pcOutput, measurement=NULL) {
    if (!is(pcOutput, "pcOutput")) {
        stop("Main model is not output from the panelcoder function", call.=FALSE)
    }

    if (!is.null(measurement)) {
        if (!is(measurement, "pcOutput")) {
            stop("Measurement model is not output from the panelcoder function", call.=FALSE)
        }
    }
    
    if (!is.null(measurement)) {
        if (!identical(pcOutput[[2]], measurement[[2]])) {
            stop("Measurement model and comparison model do not use the same data", call.=FALSE)
        }
    }

    
    if (pcOutput[[2]]$gen$yVar == TRUE) {
        maxInd <- max(pcOutput[[2]]$y$indicators,
                      pcOutput[[2]]$x$indicators)
    } else {
        maxInd <- pcOutput[[2]]$x$indicators
    }

    if (pcOutput[[2]]$gen$yVar == TRUE) {
        varNames <- c(pcOutput[[2]]$x$name, pcOutput[[2]]$y$name)
    } else {
        varNames <- pcOutput[[2]]$x$name
    }
    
    if (maxInd > 1 & is.null(measurement)) {
        warning("Only implied stability coefficients are shown for models with multiple indicators, unless a measurement model is supplied.")
        if (pcOutput[[2]]$gen$yVar == TRUE) {
            corMat <- pcOutput[[5]][1:2]
        } else {
            corMat <- pcOutput[[5]][1]
        }
        names(corMat) <- paste("Implied", varNames, sep = " ")
    } else if (maxInd > 1) {
        if (pcOutput[[2]]$gen$yVar == TRUE) {
            corMat <- pcOutput[[5]][1:2]
            mCorMat <- measurement[[5]][1:2]
        } else {
            corMat <- pcOutput[[5]][1]
            mCorMat <- measurement[[5]][1]
        }
        corMat <- cbind(corMat,mCorMat)
        names(corMat) <- c(paste("Implied",
                               varNames,
                               sep = " "),
                           paste("Observed",
                               varNames,
                               sep = " "))
    } else {
        corMat <- pcOutput[[5]]
        names(corMat) <- c(paste("Implied",
                               varNames,
                               sep = " "),
                           paste("Observed",
                               varNames,
                               sep = " "))
    }
    plotCors(corMat)
}

#' Prints all estimates 
#'
#' `panelEstimates()` takes a panelcoder output object as an argument and then
#' prints all estimates from lavaan or mplus. 
#'
#' @importFrom methods is
#' @param pcOutput An object created from running the `panelcoder()` command.
#' @export
panelEstimates <- function(pcOutput) {
    if (!is(pcOutput, "pcOutput")) {
        stop("This is not output from the panelcoder function", call.=FALSE)
    }

    output <- pcOutput[[4]]
    if (is(output, "lavaan")) {
        print(lavaan::parameterEstimates(output))
    } else {
        print(output$results$parameters)
    }
}

#' Prints formatted version of lavaan or mplus code
#'
#' `modelCode()` takes a panelcoder output object as an argument and then
#' prints corresponding lavaan or mplus code. 
#'
#' @importFrom methods is
#' @param pcOutput An object created from running the `panelcoder()` command.
#' @export
modelCode <- function(pcOutput) {
    if (!is(pcOutput, "pcOutput")) {
        stop("This is not output from the panelcoder function", call.=FALSE)
    }

    cat(pcOutput[[6]])
}



################################
#####CREATE PARCEL FUNCTION#####
################################

#' Creates parcels from items
#'
#' `parcel()` takes a set of items from a single scale and creates parcels based
#' on the user specification. The function conducts a factor analysis of items
#' at each wave and then calculates the average factor loading for each item.
#' The item with the strongest loading is assigned to the first parcel, the
#' item with the second strongest loading is assigned to the second parcel,
#' and so on until all items are assigned. 
#'
#'
#' @param data Dataframe with just the items to be parcelled.
#' @param items Number of items for the scale.
#' @param waves Number of waves.
#' @param parcels Number of parcels to create. 
#'
#' @export
parcel <- function(data,items,waves,parcels) {

    split_integer <- function(number, parts) {
        quotient <- floor(number / parts)
        remainder <- number %% parts
        
        result <- rep(quotient, times = parts)
        
        if (remainder > 0) {
            result[1:remainder] <- result[1:remainder] + 1
        }
        
        return(result)
    }

    pattern <- split_integer(items, parcels)
    
    ## Create Temporary Matrix
    loadingSum <- matrix(0,items,1)

    ## Get Average Loading
    for (i in 1:waves) {
        start=1+(i-1)*items
        finish=i*items
        loadingSum <- loadingSum+stats::factanal(na.omit(data[,start:finish]),
                                                 1)$loadings[,1]
    }
    loadingAvg <- loadingSum/waves

    ## Rank Items by Loading
    itemList <- data.frame(cbind(loadingAvg,matrix(1:items,items,1)))
    itemList <- itemList[order(itemList$X1, decreasing = TRUE),]

    create_matrix <- function(rows, cols) {
        ## Initialize an empty matrix
        mat <- matrix(0, nrow = rows, ncol = cols)
        
        ## Fill the matrix with sequential values
        counter <- 1
        for (i in 1:rows) {
            for (j in 1:cols) {
                mat[i, j] <- counter
                counter <- counter + 1
            }
        }
        
        ## Reverse the even rows in the matrix
        mat[seq(2, rows, by = 2), ] <- mat[seq(2, rows, by = 2), ncol(mat):1]
        
        return(mat)
    }

    temp <- create_matrix(max(pattern), length(pattern))

    for (row in 1:nrow(temp)) {
        for (col in 1:ncol(temp)) {
            temp[row,col] <- itemList[temp[[row,col]],2]
        }
    }

    itemMatrix <- as.data.frame(temp)
    names(itemMatrix) <- paste("parcel", 1:parcels, sep="_")
    print(itemMatrix)

    
    ##Create Parcels
    for (i in 1:waves){
        for (j in 1:parcels) {
            ##print(j)
            parcel.idx <- temp[, j][temp[, j] %in% 1:items] + ((i - 1) * items)
            ##print(parcel.idx)
            if (length(parcel.idx) > 1) {
                parcel <- apply(data[, parcel.idx], 1, FUN = mean, na.rm = TRUE)
            } else {
                parcel <- data[, parcel.idx]
            }
            index <- paste0("W", i, "P", j)
            data[,paste0("parcel",index)] <- parcel
        }
    }

    return(data)

}



################################
#####END PARCELS FUNCTION#######
################################


################################################################################
## Function to generate STARTS/CLPM/RI-CLPM Data
################################################################################
##
## Credits:
##
## Some of the initial guidance for simulating the data was taken from here:
## https://bookdown.org/marklhc/notes/simulation-example-on-structural-equation-modeling-sem.html
##
## Stationarity constraints were based on those in the STARTS model here:
## https://github.com/alexanderrobitzsch/STARTS


#' Generate STARTS Data
#'
#' `gen_starts()` generates simulated data based on the STARTS model and its
#' variants. 
#'
#' @param n Numeric value the number of lines of data to generate. Defaults to
#'   500.
#' @param nwaves Numeric value specifying the number of waves to generate.
#'   Defaults to 10.
#' @param ri_x Numeric value specifying the variance for the random intercept
#'   for X. Defaults to 1.
#' @param ri_y Numeric value specifying the variance for the random intercept
#'   for Y. Defaults to 1
#' @param cor_i Numeric value specifying the correlation between the random
#'   intercepts. Defaults to .5. 
#' @param x Numeric value specifying the variance for the autoregressive
#'   component for X. Defaults to 1.
#' @param y Numeric value specifying the variance for the autoregressive
#'   component for Y. Defaults to 1.
#' @param stab_x Numeric value specifying the stability of the autoregressive
#'   process for X. Defaults to .5.
#' @param stab_y Numeric value specifying the stability of the autoregressive
#'   process for Y. Defaults to .5.
#' @param yx Numeric value specifying the cross-lagged path predicting Y from X.
#'   Defaults to .4.
#' @param xy Numeric value specifying the cross-lagged path predicting X from Y.
#'   Defaults to .2.
#' @param cor_xy Numeric value specifying the correlation between the initial
#'   autoregressive components of X and Y. Defaults to .5.
#' @param xr Numeric value specifying the variance of the "state" or measurement
#'   error component for X. Defaults to 0.
#' @param yr Numeric value specifying the variance of the "state" or measurement
#'   error component for Y. Defaults to 0.
#' @returns Dataframe with simulated data.
#' @export
gen_starts <- function(n=500,      # N to generate
                       nwaves=10,   # Number of waves
                       ri_x=1,     # Random intercept variance for X
                       ri_y=1,     # Random intercept variance for Y
                       cor_i=.5,   # Correlation between intercepts (as correlation)
                       x=1,        # AR variance for X
                       y=1,        # AR variance for Y
                       stab_x=.5,  # Stability of X
                       stab_y=.5,  # Stability of Y
                       yx=.4,      # Cross lag (Y regressed on X)
                       xy=.2,      # Cross lag (X regressed on Y)
                       cor_xy=.5,  # Correlation between X and Y (as correlation)
                       xr=0,       # Measurement error for X
                       yr=0        # Measurement error for Y
                       ) {

    ## Transform correlations into covariances for matrices
    cor_i <- cor_i * (sqrt(ri_x) * sqrt(ri_y))
    cor_xy <- cor_xy * (sqrt(x) * sqrt(y))

    ## Stationarity Constraints
    ifelse(x==0, wxr <- 0, wxr <- (1-stab_x^2)*x - 2*stab_x*xy*cor_xy - xy^2*y)
    ifelse(y==0, wyr <- 0, wyr <- (1-stab_y^2)*y - 2*stab_y*yx*cor_xy - yx^2*x)
    ifelse(x==0 | y==0,
           cor_xyr <- 0,
           cor_xyr <- (1-stab_x*stab_y-xy*yx)*cor_xy - stab_x*yx*x - stab_y*xy*y)
    ## ifelse(x == 0 | y == 0,
    ##     cor_xyr <- 0,
    ##     cor_xyr <- cor_xy - cor_xy * (stab_x * stab_y) -
    ##         cor_xy * (xy * yx) -
    ##         (stab_x * yx * wxr) -
    ##         (stab_y * xy * wyr)
    ## )

    ## Initialize Matrices
    lambda <- matrix(0, nrow = 2 * nwaves, ncol = 2 + 2 * nwaves,
                     dimnames = list(c(paste("x",1:nwaves, sep = "_"),
                                       paste("y",1:nwaves, sep = "_")),
                                     c("ri_x", "ri_y",
                                       paste("x",1:nwaves, sep = "_"),
                                       paste("y",1:nwaves, sep = "_"))))
    theta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                    dimnames= list(c(paste("x", 1:nwaves, sep = "_"),
                                    paste("y", 1:nwaves, sep = "_")),
                                    c(paste("x", 1:nwaves, sep = "_"),
                                    paste("y", 1:nwaves, sep = "_"))))
    psi <- matrix(0, nrow = 2 + 2 * nwaves, ncol = 2 + 2 * nwaves,
                  dimnames = list(c("ri_x", "ri_y", paste("x",1:nwaves, sep = "_"),
                                    paste("y", 1:nwaves, sep = "_")),
                                  c("ri_x", "ri_y", paste("x",1:nwaves, sep = "_"),
                                    paste("y", 1:nwaves, sep = "_"))))
    beta <- matrix(0, nrow = 2 + 2 * nwaves, ncol = 2 + 2 * nwaves,
                   dimnames = list(c("ri_x", "ri_y", paste("x",1:nwaves, sep = "_"),
                                     paste("y",1:nwaves, sep = "_")),
                                   c("ri_x", "ri_y", paste("x",1:nwaves, sep = "_"),
                                     paste("y",1:nwaves, sep = "_"))))
    ##
    ## Fill in Matrices
    ## lambda
    lambda[1:nwaves, 1] <- 1 ## X loadings
    lambda[(nwaves+1):(2*nwaves), 2] <- 1  ## Y loadings
    for (i in 1:(2*nwaves)) {
        lrow <- i
        lcol <- i + 2
        lambda[lrow, lcol] <- 1
    }
    ## theta
    theta[1:nwaves, 1:nwaves] <- diag(xr, nrow = nwaves)
    theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)] <- diag(yr, nrow = nwaves)
    ## psi
    psi[1:2,1:2] <- c(ri_x, cor_i, cor_i, ri_y)
    diag(psi)[3:(2*nwaves+2)] <- c(x, rep(wxr, nwaves-1), y, rep(wyr, nwaves-1))
    psi[(nwaves+3), 3] <- cor_xy
    psi[3, (nwaves+3)] <- cor_xy
    for (i in 2:nwaves) {
        prow <- i + nwaves + 2
        pcol <- i + 2
        psi[prow, pcol] <- cor_xyr
        psi[pcol, prow] <- cor_xyr
    }
    ## beta
    for (i in 1:(nwaves-1)) {
        ## x stabilities
        xsrow <- i+3
        xscol <- i+2
        beta[xsrow, xscol] <- stab_x
        ## y stabilities
        ysrow <- i+3+nwaves
        yscol <- i+2+nwaves
        beta[ysrow, yscol] <- stab_y
    }
    for (i in 1:(nwaves-1)) {
        ## y~x cross-lagged
        ycrow <- i+3+nwaves
        yccol <- i+2
        beta[ycrow, yccol] <- yx
        ## x~y cross-lagged
        xcrow <- i+3
        xccol <- i+2+(nwaves)
        beta[xcrow, xccol] <- xy
    }
    ## Remove rows from matrices before generating data if no variance
    ## adjust psi matrix

    toDelete <- c()
    if(ri_x==0) toDelete <- c(1)
    if(ri_y==0) toDelete <- c(toDelete, 2)
    if(x==0) toDelete <- c(toDelete, c(3:(2+nwaves)))
    if(y==0) toDelete <- c(toDelete, c((3+nwaves):(3+(2*nwaves))))
    ## Delete rows and columns
    if(!is.null(toDelete)) psi <- psi[-toDelete, -toDelete]

    ## adjust beta matrix
    ## Delete rows and columns
    if(!is.null(toDelete)) beta <- beta[-toDelete, -toDelete]
        

    ## adjust lambda matrix
    ## Delete rows and columns
    if(!is.null(toDelete)) lambda <- lambda[,-toDelete]
    
    diag_length <- (x>0)*nwaves + (y>0)*nwaves + sum(ri_x>0, ri_y>0) ## Dimensions of identity matrix
    ## Generate latent factor scores
    eta <- mnormt::rmnorm(n, varcov = (solve(diag(diag_length)-beta) %*%
                               psi %*% t(solve(diag(diag_length)-beta))))
    ## Generate residuals (all zero in ri-clpm)
    ifelse(xr==0,
           ex <- matrix(0, nrow = n, ncol = nwaves),
           ex <- mnormt::rmnorm(n, varcov = theta[1:nwaves,1:nwaves]))
    ifelse(yr==0,
           ey <- matrix(0, nrow = n, ncol = nwaves),
           ey <- mnormt::rmnorm(n, varcov = theta[(nwaves+1):(2*nwaves),(nwaves+1):(2*nwaves)]))
    e <- cbind(ex,ey)
    ## Compute observed scores
    obs <- tcrossprod(eta, lambda) + e
    ## Make it a dataframe
    data <- as.data.frame(obs)
    return(data)
}

#' Add Indicators
#'
#' `addIndicators()` takes data generated using `gen_starts()` and adds
#' indicators for each measure at each wave. The amount of unique variance
#' in the indicators can be specified (though this has to be the same for all
#' indicators.
#'
#' @param df A dataframe, usually generated by `gen_starts()`. Variables should
#'   be names X1 to Xw and Y1 to Yw, where w is the last wave.
#' @param var Numeric value representing the residual variance for the
#'   indicators.
#' @param indicators Numeric value representing the number of indicators to add.
#' @returns Dataframe with the original variables and the new indicators.
#' @param labelType Either "numbers" (the default) or "letters"
#' @param sd Numeric value indicating the standard deviation of the residual
#'   variance that is added.
#' @export
addIndicators <- function(df, var, indicators, labelType="numbers", sd = 1) {
    var <- rlang::sym(var)
    if (labelType == "letters") {
        labels <- letters[1:indicators]
    } else if (labelType == "numbers") {
        labels <- paste0("_", 1:indicators)
    } else {
        stop("Label Type not recognized")
    }

    for (i in 1:indicators) {
        label <- labels[i]
        var <- rlang::enquo(var)
        prefix <- rlang::as_label(var)
        df <- df %>%
            dplyr::rowwise() %>%
            dplyr::mutate("{ prefix }{label}" := !!var + rnorm(1, 0, sd))
    }
    return(df)
}


################################################################################
## Functions to check for errors in mplus output
################################################################################

check_nonpos <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
  warning_text <- paste(
    "THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES MAY NOT BE",
    "TRUSTWORTHY FOR SOME PARAMETERS DUE TO A NON-POSITIVE DEFINITE",
    "FIRST-ORDER DERIVATIVE PRODUCT MATRIX.  THIS MAY BE DUE TO THE STARTING",
    "VALUES BUT MAY ALSO BE AN INDICATION OF MODEL NONIDENTIFICATION.",
    sep = " "
  )
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}


check_se <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
  warning_text <- paste(
    "THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES COULD NOT BE",
    "COMPUTED.  THE MODEL MAY NOT BE IDENTIFIED.  CHECK YOUR MODEL.",
    sep = " "
  )
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}

check_convergence <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
    warning_text <- paste(
        "NO CONVERGENCE.  NUMBER OF ITERATIONS EXCEEDED.",
        sep = " "
  )
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}

check_sparse <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
    warning_text <- paste(
        "THE MISSING DATA EM ALGORITHM FOR THE H1 MODEL",
        "HAS NOT CONVERGED WITH RESPECT TO THE PARAMETER",
        "ESTIMATES.  THIS MAY BE DUE TO SPARSE DATA",
        "LEADING TO A SINGULAR COVARIANCE MATRIX ESTIMATE.",
        sep = " "
  )
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}

check_coverage <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
    warning_text <- "THE COVARIANCE COVERAGE FALLS BELOW THE SPECIFIED LIMIT."
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}



#' Get estimates for lags from lavaan output
#' @param fitObject Lavaan object
#' @noRd
.summarizeLavaanLags <- function(info,
                                 fitObject,
                                 lags) {
    ## Calculate values for summary
    est <- lavaan::parameterEstimates(fitObject)
    est.std <- lavaan::standardizedSolution(fitObject)

    x_name <- info$x$name
    y_name <- info$y$name
    waves <- info$gen$maxWaves

    lagParameters <- data.frame(
        lhs = c(
            rep(paste("a", x_name, waves, sep = "_"), lags),
            rep(paste("a", y_name, waves, sep = "_"), lags),
            rep(paste("a", y_name, waves, sep = "_"), lags),
            rep(paste("a", x_name, waves, sep = "_"), lags)
        ),
        op = rep("~", 4 * lags),
        rhs = c(
            paste("a", x_name, (waves - lags):(waves - 1), sep = "_"),
            paste("a", y_name, (waves - lags):(waves - 1), sep = "_"),
            paste("a", x_name, (waves - lags):(waves - 1), sep = "_"),
            paste("a", y_name, (waves - lags):(waves - 1), sep = "_")
        )
    )

    ## Make keys to match parameters
    keyU <- with(est, paste(lhs, op, rhs, sep = "."))
    keyS <- with(est.std, paste(lhs, op, rhs, sep = "."))
    keyM <- with(lagParameters, paste(lhs, op, rhs, sep = "."))
    rows_to_keep <- match(keyM, keyU)
    rows_to_keep_s <- match(keyM, keyS)

    lagResults <- est[rows_to_keep, c("lhs", "op", "rhs", "est", "se", "pvalue")]
    lagResults.std <- est.std[rows_to_keep, "est.std"]
    lagResults <- cbind(lagResults, lagResults.std)

    lagResults$Parameter <- c(
        paste(rep(x_name, lags), "AR Lag", lags:1, sep = " "),
        paste(rep(y_name, lags), "AR Lag", lags:1, sep = " "),
        paste(rep(y_name, lags), "CL Lag", lags:1, sep = " "),
        paste(rep(x_name, lags), "CL Lag", lags:1, sep = " ")
    )

    lagResults <- lagResults[c(8, 4, 5, 6, 7)]
    names(lagResults) <- c("Parameter", "est", "se", "pvalue", "est.std")

    return(lagResults)
}


#' Summarize info from mplus object
#' @param info information object generated from `getInfo()`
#' @param fitObject mplus object with output
#' @noRd
.summarizeMplusLags <- function(info,
                                fitObject,
                                lags) {
    ## Extract estimates and fit info
    est <- fitObject$results$parameters$unstandardized
    est.std <- fitObject$results$parameters$stdyx.standardized

    x_name <- info$x$name
    x_name_up <- toupper(x_name) ## required for mplus
    y_name <- info$y$name
    y_name_up <- toupper(y_name) ## required for mplus
    waves <- info$gen$maxWaves

    lagParameters <- data.frame(
        paramHeader = c(
            rep(paste(paste("A", x_name_up, waves, sep = "_"), "ON", sep = "."), lags),
            rep(paste(paste("A", y_name_up, waves, sep = "_"), "ON", sep = "."), lags),
            rep(paste(paste("A", y_name_up, waves, sep = "_"), "ON", sep = "."), lags),
            rep(paste(paste("A", x_name_up, waves, sep = "_"), "ON", sep = "."), lags)
        ),
        param = c(
            paste("A", x_name_up, (waves - lags):(waves - 1), sep = "_"),
            paste("A", y_name_up, (waves - lags):(waves - 1), sep = "_"),
            paste("A", x_name_up, (waves - lags):(waves - 1), sep = "_"),
            paste("A", y_name_up, (waves - lags):(waves - 1), sep = "_")
        )
    )

    ## Make keys to match parameters
    keyU <- with(est, paste(paramHeader, param, sep = "."))
    keyS <- with(est.std, paste(paramHeader, param, sep = "."))
    keyM <- with(lagParameters, paste(paramHeader, param, sep = "."))
    rows_to_keep <- match(keyM, keyU)
    rows_to_keep_s <- match(keyM, keyS)

    lagResults <- est[rows_to_keep, c("paramHeader", "param", "est", "se", "pval")]
    lagResults.std <- est.std[rows_to_keep, "est"]
    lagResults <- cbind(lagResults, lagResults.std)

    lagResults$Parameter <- c(
        paste(rep(x_name, lags), "AR Lag", lags:1, sep = " "),
        paste(rep(y_name, lags), "AR Lag", lags:1, sep = " "),
        paste(rep(y_name, lags), "CL Lag", lags:1, sep = " "),
        paste(rep(x_name, lags), "CL Lag", lags:1, sep = " ")
    )

    lagResults <- lagResults[c(7, 3, 4, 5, 6)]
    names(lagResults) <- c("Parameter", "est", "se", "pvalue", "est.std")

    return(lagResults)

}


#' Print lags
#' @param program Which program used
#' @param info Basic info about model 
#' @param fitObject Lavaan object
#' @noRd
getLags <- function(pcOutput) {
    program <- pcOutput[[1]]$program
    if (program == "Mplus") {
        lags <- .summarizeMplusLags(
            info = pcOutput[[2]],
            fitObject = pcOutput[[4]],
            lags = pcOutput[[1]]$lags
        )
    } else if (program == "Lavaan") {
        lags <- .summarizeLavaanLags(
            info = pcOutput[[2]],
            fitObject = pcOutput[[4]],
            lags = pcOutput[[1]]$lags
        )
    }
    return(lags)
}


    
