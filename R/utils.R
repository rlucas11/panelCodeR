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
    if (info$gen$yVar == TRUE) {
        trait.y <- est[which(est$paramHeader == "Variances" &
                             est$param == t.y), "est"]
    } else {
        trait.y <- NA
    }
    ar.x <- est[which(est$paramHeader == "Variances" &
                         est$param == a.x), "est"]
    if (info$gen$yVar == TRUE) {
        ar.y <- est[which(est$paramHeader == "Variances" &
                          est$param == a.y), "est"]
    } else {
        ar.y <- NA
    }
    state.x <- est[which(est$paramHeader == "Variances" &
                         est$param == s.x), "est"]
    if (info$gen$yVar == TRUE) {
        state.y <- est[which(est$paramHeader == "Variances" &
                             est$param == s.y), "est"]
    } else {
        state.y <- NA
    }
    trait.x.p <- trait.x/(sum(trait.x, ar.x, state.x))
    trait.y.p <- trait.y/(sum(trait.y, ar.y, state.y))
    ar.x.p <- ar.x/(sum(trait.x, ar.x, state.x))
    ar.y.p <- ar.y/(sum(trait.y, ar.y, state.y))
    state.x.p <- state.x/(sum(trait.x, ar.x, state.x))
    state.y.p <- state.y/(sum(trait.y, ar.y, state.y))
    
    ## Correlations
    if (info$gen$yVar == TRUE) {
        trait.cor <- est.std[which(est.std$paramHeader == paste0(t.x, ".WITH") &
                                   est.std$param == t.y), "est"]
    } else {
        trait.cor <- ""
    }
    if (info$gen$yVar == TRUE) {
        ar.cor <- est.std[which(est.std$paramHeader == paste0(a.x, ".WITH") &
                                est.std$param == a.y), "est"]
    } else {
        ar.cor <- ""
    }
    if (info$gen$yVar == TRUE) {
        state.cor <- est.std[which(est.std$paramHeader == paste0(s.x, ".WITH") &
                                   est.std$param == s.y), "est"]
    } else {
        state.cor <- ""
    }
    
    ## Stability
    x.stab <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                            est.std$param == a.x), "est"]
    if (info$gen$yVar == TRUE) {
        y.stab <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                                est.std$param == a.y), "est"]
    } else {
        y.stab <- ""
    }
    
    ## Cross-lags
    if (info$gen$yVar == TRUE) {
        yOnX <- est.std[which(est.std$paramHeader == paste0(a.y2, ".ON") &
                              est.std$param == a.x), "est"]
        xOnY <- est.std[which(est.std$paramHeader == paste0(a.x2, ".ON") &
                              est.std$param == a.y), "est"]
    } else {
        yOnX <- ""
        xOnY <- ""
    }
    
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
        vars <- "x"
    }
    
    ## Get implied correlations
    if (program == "mplus") {
        impliedCors <- getMplusImpliedCors(fitObject$results$tech4$latCorEst,
                                       varNames,
                                       info$gen$maxWaves)
    } else {
        impliedCors <- getLavaanImpliedCors(fitObject,
                                            varNames,
                                            info$gen$maxWaves)
    }

    ## Get observed correlations

    if (length(vars) == 1) {
        observedCors <- summarizeR(getObservedCors(df, info, "x"))
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
        rep(c("Implied", "Observed"), each = length(vars)),
        rep(vars, 2))
    cors <- cors %>%
        dplyr::mutate(lag=row_number()) %>%
        tidyr::pivot_longer(!lag) %>%
        tidyr::separate_wider_delim(name, " ", names = c("Source", "Variable"))
    minCor <- min(cors$value)

    ggplot2::ggplot(aes(x=lag, y=value, linetype=Source, color=Variable), data=cors) +
        geom_line() +
        geom_point() +
        ylim(min(minCor,0), 1)
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
#' @param pattern Vector specifying the number of parcels and the number of
#'   items per parcel. For example, to get three parcels with 3 items for the
#'   first parcel, 3 items for the second, and 4 items for the third, specify
#'   `pattern = c(3, 3, 4)`
#'
#' @export
parcel <- function(data,items,waves,pattern) {

#Create Temporary Matrix
loadingSum <- matrix(0,items,1)

#Get Average Loading
for (i in 1:waves) {
  start=1+(i-1)*items
  finish=i*items
  loadingSum <- loadingSum+stats::factanal(na.omit(data[,start:finish]),
                                           1)$loadings[,1]
}
loadingAvg <- loadingSum/waves

#Rank Items by Loading
itemList <- data.frame(cbind(loadingAvg,matrix(1:items,items,1)))
itemList <- itemList[order(itemList$X1),]

#Create Pattern Matrix for Creating Parcels
temp <- matrix(nrow=max(pattern),ncol=length(pattern))
for (i in 1:max(pattern)) {
  if (gtools::odd(i)) {
    temp[i,] <- ((1+((i-1)*(length(pattern)))):(3+((i-1)*(length(pattern)))))
  }
  else {
    temp[i,] <- ((length(pattern)*i):(length(pattern)*i-2))
  }
}

#Create Parcels
for (i in 1:waves){
  if (gtools::odd(max(pattern))) newPattern <- pattern else {
    newPattern <- pattern
    newPattern[1] <- pattern[3]
    newPattern[2] <- pattern[2]
    newPattern[3] <- pattern[1]
  }
  for (j in 1:length(newPattern)) {
    parcel=0
    maxrows <- newPattern[j]  
    for (k in 1:maxrows) {
      parcel <- parcel+data[,((i-1)*items+itemList[temp[k,j],2])]
    }
    indexw <- paste("W",i,sep="")
    indexp <- paste("P",j,sep="")
    index <- paste(indexw,indexp,sep="")
    data[,paste("parcel",index,sep="")] <- parcel/newPattern[j]
  }
}

return(data[,(items*waves+1):(items*waves+length(pattern)*waves)])

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
    ## ifelse(x==0 | y==0,
    ##        cor_xyr <- 0,
    ##        cor_xyr <- (1-stab_x*stab_y-xy*yx)*cor_xy - stab_x*yx*x - stab_y*xy*y)
    ifelse(x == 0 | y == 0,
        cor_xyr <- 0,
        cor_xyr <- cor_xy - cor_xy * (stab_x * stab_y) -
            cor_xy * (xy * yx) -
            (stab_x * yx * wxr) -
            (stab_y * xy * wyr)
    )

    ## Initialize Matrices
    lambda <- matrix(0, nrow = 2 * nwaves, ncol = 2 + 2 * nwaves,
                     dimnames = list(c(paste0("x",1:nwaves),
                                       paste0("y",1:nwaves)),
                                     c("ri_x", "ri_y", paste0("x",1:nwaves),
                                       paste0("y",1:nwaves))))
    theta <- matrix(0, nrow = 2 * nwaves, ncol = 2 * nwaves,
                    dimnames= list(c(paste0("x", 1:nwaves),
                                    paste0("y", 1:nwaves)),
                                    c(paste0("x", 1:nwaves),
                                    paste0("y", 1:nwaves))))
    psi <- matrix(0, nrow = 2 + 2 * nwaves, ncol = 2 + 2 * nwaves,
                  dimnames = list(c("ri_x", "ri_y", paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves)),
                                  c("ri_x", "ri_y", paste0("x",1:nwaves),
                                    paste0("y", 1:nwaves))))
    beta <- matrix(0, nrow = 2 + 2 * nwaves, ncol = 2 + 2 * nwaves,
                   dimnames = list(c("ri_x", "ri_y", paste0("x",1:nwaves),
                                     paste0("y",1:nwaves)),
                                   c("ri_x", "ri_y", paste0("x",1:nwaves),
                                     paste0("y",1:nwaves))))
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
