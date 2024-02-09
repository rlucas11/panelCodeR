
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
                waves = xWaves,
                actualWaves = xWaves,
                indicators = xIndicators
            ),
            y = NULL
        )
    }
    return(info)
}


.summarizeLavaan <- function(fitObject) {
    ## Calculate values for summary
    parEst <- lavaan::parameterEstimates(fitObject)
    stdEst <- lavaan::standardizedSolution(fitObject)
    fit <- lavaan::fitMeasures(fitObject)
    ## Variance Decomp
    trait.x <- parEst[which(parEst$label=="x_tVar"), "est"]
    trait.y <- parEst[which(parEst$label=="y_tVar"), "est"]
    ar.x <- parEst[which(parEst$label=="xvar1"), "est"]
    ar.y <- parEst[which(parEst$label=="yvar1"), "est"]
    state.x <- parEst[which(parEst$label=="sx1"), "est"]
    state.y <- parEst[which(parEst$label=="sy1"), "est"]
    trait.x.p <- trait.x/(sum(trait.x, ar.x, state.x))
    trait.y.p <- trait.y/(sum(trait.y, ar.y, state.y))
    ar.x.p <- ar.x/(sum(trait.x, ar.x, state.x))
    ar.y.p <- ar.y/(sum(trait.y, ar.y, state.y))
    state.x.p <- state.x/(sum(trait.x, ar.x, state.x))
    state.y.p <- state.y/(sum(trait.y, ar.y, state.y))
    ## Correlations
    trait.cor <- stdEst[which(stdEst$label=="cov_txty"), "est.std"]
    ar.cor <- stdEst[which(stdEst$label=="cov_ar1"), "est.std"]
    state.cor <- stdEst[which(stdEst$label=="cov_s1"), "est.std"]
    ## Stability
    x.stab <- parEst[which(parEst$label=="a2"), "est"]
    y.stab <- parEst[which(parEst$label=="b2"), "est"]
    ## Cross-lags
    yOnX <- parEst[which(parEst$label=="c2"), "est"]
    xOnY <- parEst[which(parEst$label=="d2"), "est"]
    ## Fit
    chi2 <- fit[["chisq"]]
    chi2df <- fit[["df"]]
    chi2p <- fit[["pvalue"]]
    cfi <- fit[["cfi"]]
    tli <- fit[["tli"]]
    srmr <- fit[["srmr"]]
    rmsea <- fit[["rmsea"]]
    outputList <- list(
        trait.x = trait.x.p,
        trait.y = trait.y.p,
        ar.x = ar.x.p,
        ar.y = ar.y.p,
        state.x = state.x.p,
        state.y = state.y.p,
        trait.cor = trait.cor,
        ar.cor = ar.cor,
        state.cor = state.cor,
        x.stab = x.stab,
        y.stab = y.stab,
        yOnX = yOnX,
        xOnY = xOnY,
        chi2 = chi2,
        chi2df = chi2df,
        chi2p = chi2p,
        cfi = cfi,
        tli = tli,
        srmr = srmr,
        rmsea = rmsea)
    class(outputList) <- "pcSum"
    return(outputList)
}


.summarizeMplus <- function(info, fitObject) {
    ## Extract estimates and fit info
    est <- fitObject$results$parameters$unstandardized
    est.std <- fitObject$results$parameters$stdyx.standardized
    fit <- fitObject$results$summaries

    ## Set names
    t.x  <- paste("T", info$x$name, sep = "_")
    if (info$gen$yVar == TRUE) {
        t.y <- paste("T", info$y$name, sep = "_")
    }
    a.x  <- paste("A", info$x$name, "1", sep = "_")
    a.x2  <- paste("A", info$x$name, "2", sep = "_")
    if (info$gen$yVar == TRUE) {
        a.y <- paste("A", info$y$name, "1", sep = "_")
        a.y2 <- paste("A", info$y$name, "2", sep = "_")
    }
    s.x  <- paste("S", info$x$name, "1", sep = "_")
    if (info$gen$yVar == TRUE) {
        s.y <- paste("S", info$y$name, "1", sep = "_")
    }
    
    ## Variance Decomp
    trait.x <- est[which(est$paramHeader == "Variances" &
                         est$param == t.x), "est"]
    trait.y <- est[which(est$paramHeader == "Variances" &
                         est$param == t.y), "est"]
    ar.x <- est[which(est$paramHeader == "Variances" &
                         est$param == a.x), "est"]
    ar.y <- est[which(est$paramHeader == "Variances" &
                         est$param == a.y), "est"]
    state.x <- est[which(est$paramHeader == "Variances" &
                         est$param == s.x), "est"]
    state.y <- est[which(est$paramHeader == "Variances" &
                         est$param == s.y), "est"]
    trait.x.p <- trait.x/(sum(trait.x, ar.x, state.x))
    trait.y.p <- trait.y/(sum(trait.y, ar.y, state.y))
    ar.x.p <- ar.x/(sum(trait.x, ar.x, state.x))
    ar.y.p <- ar.y/(sum(trait.y, ar.y, state.y))
    state.x.p <- state.x/(sum(trait.x, ar.x, state.x))
    state.y.p <- state.y/(sum(trait.y, ar.y, state.y))
    
    ## Correlations
    trait.cor <- est.std[which(est.std$paramHeader == paste0(t.x, ".WITH") &
                               est.std$param == t.y), "est"]
    ar.cor <- est.std[which(est.std$paramHeader == paste0(a.x, ".WITH") &
                               est.std$param == a.y), "est"]
    state.cor <- est.std[which(est.std$paramHeader == paste0(s.x, ".WITH") &
                               est.std$param == s.y), "est"]
    
    ## Stability
    x.stab <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                            est.std$param == a.x), "est"]
    y.stab <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                            est.std$param == a.y), "est"]
    
    ## Cross-lags
    yOnX <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                          est.std$param == a.x), "est"]
    xOnY <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                          est.std$param == a.y), "est"]
    
    ## Fit
    chi2 <- fit[["ChiSqM_Value"]]
    chi2df <- fit[["ChiSqM_DF"]]
    chi2p <- fit[["ChiSqM_PValue"]]
    cfi <- fit[["CFI"]]
    tli <- fit[["TLI"]]
    srmr <- fit[["SRMR"]]
    rmsea <- fit[["RMSEA_Estimate"]]
    outputList <- list(
        trait.x = trait.x.p,
        trait.y = trait.y.p,
        ar.x = ar.x.p,
        ar.y = ar.y.p,
        state.x = state.x.p,
        state.y = state.y.p,
        trait.cor = trait.cor,
        ar.cor = ar.cor,
        state.cor = state.cor,
        x.stab = x.stab,
        y.stab = y.stab,
        yOnX = yOnX,
        xOnY = xOnY,
        chi2 = chi2,
        chi2df = chi2df,
        chi2p = chi2p,
        cfi = cfi,
        tli = tli,
        srmr = srmr,
        rmsea = rmsea)
    class(outputList) <- "pcSum"
    return(outputList)
}

#' Summarizes results from `panelcoder()`
#'
#' @param object Results from `panelcoder()`.
#' @param ... Additional arguments to `summary()`.
#'
#' @export
summary.pcSum <- function(object, ...) {
    cat("Model Summary: \n")
    cat("Model: Chi2 (df = ",
        sprintf("%i", object$chi2df), ") = ",
        sprintf("%.3f", object$chi2), ", p = ",
        sprintf("%.3f", object$chi2p), "\n",
        sep="")
    cat("\n")
    cat("Fit Indices:")
    cat("\n")
    cat("CFI = ",
        sprintf("%.3f", object$cfi),
        ", TLI = ",
        sprintf("%.3f", object$tli),
        ", SRMR = ",
        sprintf("%.3f", object$srmr), "\n",
        sep="")
    cat("RMSEA = ",
        sprintf("%.3f", object$rmsea),
        "\n")
    cat("\n")
    cat("Variance Decomposition: \n")
    cat("\n")
    cat("X Variable: \n")
    cat("Trait: ",
        sprintf("%.3f", object$trait.x),
        ", Autoregressive: ",
        sprintf("%.3f", object$ar.x),
        ", State: ",
        sprintf("%.3f", object$state.x),
        "\n",
        sep="")
    cat("Y Variable: \n")
    cat("Trait: ",
        sprintf("%.3f", object$trait.y),
        ", Autoregressive: ",
        sprintf("%.3f", object$ar.y),
        ", State: ",
        sprintf("%.3f", object$state.y), "\n",
        sep="")
    cat("\n")
    cat("Stability: \n")
    cat("X: ", sprintf("%.3f", object$x.stab), "\n")
    cat("Y: ", sprintf("%.3f", object$y.stab), "\n")
    cat("\n")
    cat("Cross-Lag Paths: \n")
    cat("Y predicted from X: ", sprintf("%.3f", object$yOnX) ,"\n", sep="")
    cat("X predicted from Y: ", sprintf("%.3f", object$xOnY),"\n", sep="")
    cat("\n")
    cat("Correlations: \n")
    cat("\n")
    cat("Stable Trait: ", sprintf("%.3f", object$trait.cor), "\n")
    cat("Autoregressive Trait: ", sprintf("%.3f", object$ar.cor), "\n")
    cat("State: ", sprintf("%.3f", object$state.cor), "\n")
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


combineCors <- function(df, info, program, fitObject) {
    if (info$gen$yVar == TRUE) {
        varNames <- c(info$x$name, info$y$name)
        vars <- c("x", "y")
    } else {
        varNames <- c(info$x$name)
        vars <- c("x", "y")
    }
    
    ## Get implied correlations
    if (program == "mplus") {
        impliedCors <- getMplusImpliedCors(fitObject[[4]]$results$tech4$latCorEst,
                                       varNames,
                                       info$gen$maxWaves)
    } else {
        impliedCors <- getLavaanImpliedCors(fitObject[[4]],
                                            varNames,
                                            info$gen$maxWaves)
    }

    ## Get observed correlations

    if (length(vars) == 1) {
        observedCors <- getObservedCors(df, info, "x")
    } else {
        observedCors <- cbind(summarizeR(getObservedCors(df, info, "x")),
                              summarizeR(getObservedCors(df, info, "y")))
    }

    allCors <- cbind(impliedCors, observedCors)
    return(allCors)
}

    

plotCors <- function(cors, vars) {
    cors <- as.data.frame(cors)
    lags <- nrow(cors)
    names(cors) <- paste(
        rep(c("Implied", "Observed"), each = 2),
        rep(vars, 2))
    cors <- cors %>%
        mutate(lag=row_number()) %>%
        pivot_longer(!lag) %>%
        separate_wider_delim(name, " ", names = c("Source", "Variable"))
    minCor <- min(cors$value)

    ggplot(aes(x=lag, y=value, linetype=Source, color=Variable), data=cors) +
        geom_line() +
        geom_point() +
        ylim(min(minCor,0), 1)
}
