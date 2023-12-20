################################################################################
## Function to create and run mplus STARTS models
## Can include phantom variables for missing waves
## 
## "waves" = total number of possible waves (e.g., 20 for HILDA, 37 for SOEP)
## "actual waves" = waves with data (e.g., c(1,2,3,5,7,9))
################################################################################

################################################################################
## Build trait part of input file
################################################################################
buildTraitX <- function(waves) {
    wavesList <- 1:waves
    title <- "!!! Stable Trait X;\n"
    ## First do "BY" statements
    traitModel <- paste0(title,"traitX by lx", wavesList[1],"@1;\n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(traitModel,
                         paste0("traitX by lx", w, "@1;\n")
                         )
    }
    ## Then set variance
    traitModel <- paste0(traitModel, "\n!!! Trait Variance for X \ntraitX (tx);\n")
    return(traitModel)
}

buildTraitY <- function(waves) {
    wavesList <- 1:waves
    title <- "!!! Stable Trait Y;\n"
    ## First do "BY" statements
    traitModel <- paste0(title,"traitY by ly", wavesList[1],"@1;\n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(traitModel,
                         paste0("traitY by ly", w, "@1;\n")
                         )
    }
    ## Then set variance
    traitModel <- paste0(traitModel, "\n!!! Trait Variance for Y \ntraitY (ty);\n")
    return(traitModel)
}


################################################################################
## Build autoregressive part of input file
################################################################################

buildARX <- function(waves) {
    wavesList <- 1:waves
    title <- "!!! Autoregressive Part for X;\n"
    subtitle1 <- "\n!!! Indicator Statements;\n"
    ## First do "BY" statements linking AR to latent occasion variable
    arModel <- paste0(title, subtitle1, "arx1 by lx1@1;\n")
    for (w in wavesList[-1]) {
        arModel <- paste0(arModel, "arx", w, " by lx", w, "@1;\n")
    }
    ## Then add stability paths
    subtitle2 <- "\n!!! Regression Statements;\n"
    arModel <- paste0(arModel, subtitle2, "arx2 on arx1(a);\n")
    for (w in wavesList[-c(1:2)]) {
        arModel <- paste0(
            arModel,
            paste0("arx", w, " on arx", (w - 1), "(a", (w - 1), ");\n")
        )
    }
    ## Finally add variances (waves after the first are set to be equal)
    arModel <- paste0(
        arModel,
        "\n!!! Autoregressive Component Variance for X\narx1 (arvx);\n"
    )
    for (w in wavesList[-1]) {
        arModel <- paste0(
            arModel,
            paste0("arx", w, " (arvx", w, ");\n")
        )
    } 
    return(arModel)
}


buildARY <- function(waves) {
    wavesList <- 1:waves
    title <- "!!! Autoregressive Part for Y;\n"
    subtitle1 <- "\n!!! Indicator Statements;\n"
    ## First do "BY" statements linking AR to latent occasion variable
    arModel <- paste0(title, subtitle1, "ary1 by ly1@1;\n")
    for (w in wavesList[-1]) {
        arModel <- paste0(arModel, "ary", w, " by ly", w, "@1;\n")
    }
    ## Then add stability paths
    subtitle2 <- "\n!!! Regression Statements;\n"
    arModel <- paste0(arModel, subtitle2, "ary2 on ary1(b);\n")
    for (w in wavesList[-c(1:2)]) {
        arModel <- paste0(
            arModel,
            paste0("ary", w, " on ary", (w - 1), "(b", (w - 1), ");\n")
        )
    }
    ## Finally add variances (waves after the first are set to be equal)
    arModel <- paste0(
        arModel,
        "\n!!! Autoregressive Component Variance for Y\nary1 (arvy);\n"
    )
    for (w in wavesList[-1]) {
        arModel <- paste0(
            arModel,
            paste0("ary", w, " (arvy", w, "); \n")
        )
    }
    return(arModel)
}


################################################################################
## Setup phantom variables for missing waves (if needed; returns comment if none)
################################################################################

buildPhantomX <- function(waves, xWaves) {
    if(length(xWaves) < length(c(1:waves))) {
        title <- "!!! Phantom X Variables;\n"
        phantomX <- title
        phantomWaves <- c(1:waves)[-xWaves]
        for (w in phantomWaves) {
            phantomX <- paste0(phantomX, "lx", w, " by ;\n")
        }
        return(phantomX)
    } else {
        return("!!! No Phantom X Variables;\n")
    }
}

buildPhantomY <- function(waves, yWaves) {
    if(length(yWaves) < length(c(1:waves))) {
        title <- "!!! Phantom Y Variables;\n"
        phantomY <- title
        phantomWaves <- c(1:waves)[-yWaves]
        for (w in phantomWaves) {
            phantomY <- paste0(phantomY, "ly", w, " by ;\n")
        }
        return(phantomY)
    } else {
        return("!!! No Phantom Y Variables;\n")
    }
}



################################################################################
## Build observed part of model
################################################################################

buildObservedX <- function(xWaves, xIndicators) {
    title <- "!!! Observed Variables for X;\n"
    ## First do "BY" statements
    ## If only 1 indicator
    if (xIndicators == 1) {
        observedModel <- paste0(title, "lx", xWaves[1], " by x", xWaves[1], "@1;\n")
        for (w in xWaves[-1]) {
            observedModel <- (paste0(observedModel,
                                     "lx", w, " by x", w, "@1;\n"))
        }
        ## Then constrain variance to 0 (we have a separate state variable)
        observedModel <- paste0(observedModel,
                                "\n!!! Residual Variance Constrained to 0\nx",
                                xWaves[1], "@0;\n")
        for (w in xWaves[-1]) {
            observedModel <- paste0(observedModel,
                                    "x", w, "@0;\n")
        }
    } else {
        indLabels <- letters[1:xIndicators]
        observedModel <- "!! Multiple Indicators; \n"
        for (w in xWaves) {
            for (i in indLabels) {
                if (i == "a") {
                    observedModel <- (
                        paste0(observedModel,
                               "lx",
                               w,
                               " by x",
                               w,
                               i,
                               "@1;\n")
                    )
                } else {
                    observedModel <- (
                        paste0(observedModel,
                               "lx",
                               w,
                               " by x",
                               w,
                               i,
                               " (lx",
                               i,
                               ") ;\n")
                    )
                }
            }
        }
        ## Set residual variances for indicators
        observedModel <- paste0(observedModel,
                                "\n!!! Indicator Residual Variance \n")
        for (w in xWaves) {
            for (i in indLabels) {
                observedModel <- paste0(observedModel,
                                        "x",
                                        w,
                                        i,
                                        " (xr",
                                        i,
                                        ") ;\n")
            }
        }
    }
    return(observedModel)
}

buildObservedY <- function(yWaves, yIndicators) {
    title <- "!!! Observed Variables for Y;\n"
    ## First do "BY" statements
    ## If only 1 indicator
    if (yIndicators == 1) {
        observedModel <- paste0(title, "ly", yWaves[1], " by y", yWaves[1], "@1;\n")
        for (w in yWaves[-1]) {
            observedModel <- (paste0(observedModel,
                                     "ly", w, " by y", w, "@1;\n"))
        }
        ## Then constrain variance to 0 (we have a separate state variable)
        observedModel <- paste0(observedModel,
                                "\n!!! Residual Variance Constrained to 0\ny",
                                yWaves[1], "@0;\n")
        for (w in yWaves[-1]) {
            observedModel <- paste0(observedModel,
                                    "y", w, "@0;\n")
        }
    } else {
        indLabels <- letters[1:yIndicators]
        observedModel <- "!! Multiple Indicators; \n"
        for (w in yWaves) {
            for (i in indLabels) {
                if (i == "a") {
                    observedModel <- (
                        paste0(observedModel,
                               "ly",
                               w,
                               " by y",
                               w,
                               i,
                               "@1;\n")
                    )
                } else {
                    observedModel <- (
                        paste0(observedModel,
                               "ly",
                               w,
                               " by y",
                               w,
                               i,
                               " (ly",
                               i,
                               ") ;\n")
                    )
                }
            }
        }
        ## Set residual variances for indicators
        observedModel <- paste0(observedModel,
                                "\n!!! Indicator Residual Variance \n")
        for (w in yWaves) {
            for (i in indLabels) {
                observedModel <- paste0(observedModel,
                                        "y",
                                        w,
                                        i,
                                        " (yr",
                                        i,
                                        ") ;\n")
            }
        }
    }
    return(observedModel)
}




################################################################################
## Latent State Variance
################################################################################

buildStateX <- function(waves) {
    wavesList <- 1:waves
    title <- "\n!!! States for X;\n"
    ## First do "BY" statements
    stateModel <- paste0(title, "sx1 by lx1@1;\n")
    for (w in wavesList[-1]) {
        stateModel <- paste0(
            stateModel,
            paste0("sx", w, " by lx", w, "@1;\n")
        )
    }
    ## Then add variances
    subtitle <- "\n!!! State Variance for X;\n"
    stateModel <- paste0(stateModel, subtitle, "sx1 (sx);\n")
    for (w in wavesList[-1]) {
        stateModel <- paste0(
            stateModel,
            paste0("sx", w, " (sx", w, "); \n")
        )
    }
    return(stateModel)
}

buildStateY <- function(waves) {
    wavesList <- 1:waves
    title <- "\n!!! States for Y;\n"
    ## First do "BY" statements
    stateModel <- paste0(title, "sy1 by ly1@1;\n")
    for (w in wavesList[-1]) {
        stateModel <- paste0(stateModel,
                             paste0("sy", w, " by ly", w, "@1;\n"))
    }
    ## Then add variances
    subtitle <- "\n!!! State Variance for Y;\n"
    stateModel <- paste0(stateModel, subtitle, "sy1 (sy);\n")
    for (w in wavesList[-1]) {
        stateModel <- paste0(
            stateModel,
            paste0("sy", w, " (sy", w, "); \n")
        )
    }
    return(stateModel)
}


################################################################################
## latent occasion Variance
################################################################################

## Build latent occasion variance section
buildLatentVarX <- function(waves) {
    wavesList <- 1:waves
    ## All we need to do is set variances here, as we have defined these elsewhere
    title <- "\n!!! Constrain Latent Occasion Residuals for X to 0\n"
    latentVar <- paste0(title, "lx1@0;\n")
    for (w in wavesList[-1]) {
        latentVar <- paste0(latentVar,
                            paste0("lx", w, "@0;\n"))
    }
    return(latentVar)
}

## Build latent occasion variance section
buildLatentVarY <- function(waves) {
    wavesList <- 1:waves
    ## All we need to do is set variances here, as we have defined these elsewhere
    title <- "\n!!! Constrain Latent Occasion Residuals for Y to 0\n"
    latentVar <- paste0(title, "ly1@0;\n")
    for (w in wavesList[-1]) {
        latentVar <- paste0(latentVar,
                            paste0("ly", w, "@0;\n"))
    }
    return(latentVar)
}

################################################################################
## Cross-Lagged Paths
################################################################################

buildClYFromX <- function(wavesList) {
    title <- "\n!!! Cross-Lagged Paths \n"
    subtitle <- "\n!!! Y Predicted from X \n"
    clYonX <- paste0(title, subtitle)
    clYonX <- paste(clYonX,
        "ary2 on arx1 (c); \n",
        sep = "\n"
    )
    for (w in wavesList[-c(1:2)]) {
        clYonX <- paste0(
            clYonX, "ary", w, " on arx", (w - 1), "(c", (w - 1), "); \n"
        )
    }
    return(clYonX)
}


buildClXFromY <- function(wavesList) {
    title <- "\n!!! Cross-Lagged Paths \n"
    subtitle <- "\n!!! X Predicted from Y \n"
    clXonY <- paste0(title, subtitle)
    clXonY <- paste(clXonY,
        "arx2 on ary1 (d); \n",
        sep = "\n"
    )
    for (w in wavesList[-c(1:2)]) {
        clXonY <- paste0(
            clXonY, "arx", w, " on ary", (w - 1), "(d", (w - 1), "); \n"
        )
    }
    return(clXonY)
}

buildTraitCors <- function() {
    traitCors <- "\n!!! Stable Trait Correlations \ntraitX with traitY(cov_txty); \n"
    return(traitCors)
}

buildArCors <- function(wavesList) {
    arCors <- "!!! AR Correlations \narx1 with ary1(cov_xy); \n"
    for (w in wavesList[-1]) {
        arCors <- paste0(
            arCors,
            paste0(
                "arx",
                w,
                " with ary",
                w,
                "(cov_xyr",
                w,
                "); \n"
            )
        )
    }
    return(arCors)
}


buildStateCors <- function(wavesList) {
    stateCors <- "!!! State Correlations \n"
    stateCors <- paste(stateCors,
        "sx1 with sy1(cov_s); \n",
        sep = " \n"
    )
    for (w in wavesList[-1]) {
        stateCors <- paste0(
            stateCors,
            paste0("sx",
                   w,
                   " with sy",
                   w,
                   "(cov_s",
                   w,
                   "); \n")
        )
    }
    return(stateCors)
}


buildResidCorsX <- function(xWaves, xIndicators) {
    residCors <- "!!! Residual Correlations \n"
    indLabels <- letters[1:xIndicators]
    for (i in indLabels) {
        for (j in 1:length(xWaves)) {
            for (k in xWaves[-c(1:j)]) {
                residCors <- paste0(
                    residCors,
                    "x",
                    xWaves[j],
                    i,
                    " with x",
                    k,
                    i,
                    " (x",
                    i,
                    "lag",
                    k - xWaves[j],
                    "); \n"
                )
            }
        }
    }
    return(residCors)
}

buildResidCorsY <- function(yWaves, yIndicators) {
    residCors <- "!!! Residual Correlations \n"
    indLabels <- letters[1:yIndicators]
    for (i in indLabels) {
        for (j in 1:length(yWaves)) {
            for (k in yWaves[-c(1:j)]) {
                residCors <- paste0(
                    residCors,
                    "y",
                    yWaves[j],
                    i,
                    " with y",
                    k,
                    i,
                    " (y",
                    i,
                    "lag",
                    k - yWaves[j],
                    "); \n"
                )
            }
        }
    }
    return(residCors)
}



################################################################################
## Build final MODEL statement
################################################################################
#' Builds Mplus Code
#'
#' `buildMplus()` produces code for variations of the general Stable Trait,
#' Autoregressive Trait, State Model. Univariate and bivariate versions can
#' be specified. Reduced versions of the model can be specified by setting
#' variance components for trait, AR, or state to be `FALSE`.
#'
#' @param waves Numeric value indicating total number of waves.
#' @param XVar Logical. Include X variable.
#' @param YVar logical. Include Y variable.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param trait Logical value indicating whether to include a stable trait.
#'   Defaults to TRUE.
#' @param AR logical value indicating whether to include autoregressive trait.
#'   Defaults to TRUE.
#' @param state logical value indicating whether to include state. Defaults to
#'   TRUE.
#' @param crossLag logical value indicating whether to include cross-lagged
#'   paths. Defaults to TRUE.
#' @param stateCor logical value indicating whether to include state
#'   correlations. Defaults to TRUE.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
buildMplus <- function(waves,
                       XVar = TRUE,
                       YVar = TRUE,
                       xWaves = NULL,
                       yWaves = NULL,
                       xIndicators = 1,
                       yIndicators = 1,
                       trait = TRUE,
                       AR = TRUE,
                       state = TRUE,
                       crossLag = TRUE,
                       stateCor = FALSE,
                       stationarity = TRUE) {
    wavesList <- 1:waves
    if (is.null(xWaves)) xWaves <- wavesList
    if (is.null(yWaves)) yWaves <- wavesList
    ## Makes sure xWaves and yWaves <= total waves
    xWaves <- xWaves[xWaves <= waves]
    yWaves <- yWaves[yWaves <= waves]
    if (XVar == FALSE | YVar == FALSE) {
        crossLag <- FALSE
        stateCor <- FALSE
    }
    if (AR == FALSE) {
        crossLag <- FALSE
    }
    
    ## Start model statement
    modelStatement <- "!!! Automated Code for Mplus \n"
    
    ## Phantom variables
    if (XVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildPhantomX(waves, xWaves),
            buildObservedX(xWaves, xIndicators),
            buildLatentVarX(waves),
            buildTraitX(waves),
            buildARX(waves),
            buildStateX(waves),
            sep = " \n"
        )
    }

    if (YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildPhantomY(waves, yWaves),
            buildObservedY(yWaves, yIndicators),
            buildLatentVarY(waves),
            buildTraitY(waves),
            buildARY(waves),
            buildStateY(waves),
            sep = " \n"
        )
    }
        
    
    ## Add Cross Lags
    if (crossLag == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildClYFromX(wavesList),
            buildClXFromY(wavesList),
            sep = " \n"
        )
    }

    ## Add Trait Cors
    if (XVar == TRUE & YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildTraitCors(),
            sep = " \n"
        )
    }

    ## Add AR Cors
    if (XVar == TRUE & YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildArCors(wavesList),
            sep = " \n"
        )
    }

    ## Add State Cors
    if (XVar == TRUE & YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildStateCors(wavesList),
            sep = " \n"
        )
    }

    ## Add Resid Cors
    if (XVar == TRUE & xIndicators > 1) {
        modelStatement <- paste(
            modelStatement,
            buildResidCorsX(xWaves, xIndicators),
            sep = " \n"
        )
    }

    if (YVar == TRUE & yIndicators > 1) {
        modelStatement <- paste(
            modelStatement,
            buildResidCorsY(yWaves, yIndicators),
            sep = " \n"
        )
    }

    ## Return model
    return(modelStatement)
}


## Function to specify constraints

#' Build Constraints
#'
#' Builds constraint statement for Mplus input file based on specified options.
#' 
#' @param waves Numeric value indicating total number of waves.
#' @param XVar Logical. Include X variable.
#' @param YVar logical. Include Y variable.
#' @param trait Logical value indicating whether to include a stable trait.
#'   Defaults to TRUE. 
#' @param AR logical value indicating whether to include autoregressive trait.
#'   Defaults to TRUE.
#' @param state logical value indicating whether to include state. Defaults to
#'   TRUE.
#' @param crossLag logical value indicating whether to include cross-lagged
#'   paths. Defaults to TRUE.
#' @param stateCor logical value indicating whether to include state
#'   correlations. Defaults to TRUE.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @returns Character vector representing the Mplus code for the constraints
#'   statement.
#' @export
buildConstraints <- function(waves,
                             XVar = TRUE,
                             YVar = TRUE,
                             trait = TRUE,
                             AR = TRUE,
                             state = TRUE,
                             crossLag = TRUE,                             
                             stateCor = TRUE,
                             stationarity = TRUE,
                             constrainCor = stationarity,
                             limits = TRUE) {
    wavesList <- 1:waves
    if (XVar == FALSE | YVar == FALSE | AR == FALSE) {
        crossLag <- FALSE
    }
    ## First make constraints for different models
    ## Variances and Covariances
    modelConstraints <- " \n!! Use constraints to specify Models; \n"
    ## Set unused parameters to 0
    if (XVar == TRUE & YVar == TRUE & crossLag == FALSE) {
        modelConstraints <- paste(modelConstraints,
            "NEW(c); ",
            "c = 0; ",
            "NEW(d); ",
            "d = 0; ",
            sep = " \n"
        )
    }
    if (XVar == TRUE & YVar == FALSE) {
        modelConstraints <- paste(modelConstraints,
            "NEW(d); ",
            "d = 0; ",
            "NEW(cov_xy); ",
            "cov_xy = 0; ",
            sep = " \n"
        )
    }
    if (XVar == FALSE & YVar == TRUE) {
        modelConstraints <- paste(modelConstraints,
            "NEW(c); ",
            "c = 0; ",
            "NEW(cov_xy); ",
            "cov_xy = 0; ",
            sep = " \n"
        )
    }
    if (XVar == TRUE & AR == FALSE) {
        modelConstraints <- paste(modelConstraints,
##            "NEW(a); ",
            "a = 0; ",
##            "NEW(arvx); ",
            "arvx = 0; ",
            sep = " \n"
        )
    }
    if (YVar == TRUE & AR == FALSE) {
        modelConstraints <- paste(modelConstraints,
##            "NEW(b); ",
            "b = 0; ",
##            "NEW(arvy); ",
            "arvy = 0; ",
            sep = " \n"
        )
    }
    if ((XVar == TRUE & AR == TRUE & stationarity == TRUE) |
        (XVar == TRUE & AR == FALSE)) {
        for (w in wavesList[-c(1, waves)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "a",
                w,
                " = a; \n"
            )
        }
    }
    if ((YVar == TRUE & AR == TRUE & stationarity == TRUE) |
        (YVar == TRUE & AR == FALSE)) {
        for (w in wavesList[-c(1, waves)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "b",
                w,
                " = b; \n"
            )
        }
    }
    if (crossLag == TRUE & stationarity == TRUE) {
        for (w in wavesList[-c(1, waves)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "c",
                w,
                " = c; \n",
                "d",
                w,
                " = d; \n"
            )
        }
    }
    ## Stable Traits
    if (trait == FALSE) {
        if (XVar == TRUE) {
            modelConstraints <- paste(modelConstraints,
                "tx = 0;",
                sep = " \n"
            )
        }
        if (YVar == TRUE) {
            modelConstraints <- paste(modelConstraints,
                "ty = 0;",
                sep = " \n"
            )
        }
        if (XVar == TRUE & YVar == TRUE) {
            modelConstraints <- paste(modelConstraints,
                "cov_txty = 0; \n",
                sep = " \n"
            )
        }
    }
    
    if (state == FALSE) {
        if (XVar == TRUE) {
            modelConstraints <- paste(modelConstraints,
                "sx = 0;",
                sep = " \n"
                )
        }
        if (YVar == TRUE) {
            modelConstraints <- paste(modelConstraints,
                "sy = 0;",
                sep = " \n"
            )
        }
        if (XVar == TRUE & YVar == TRUE) {
            modelConstraints <- paste(modelConstraints,
                "cov_s = 0;",
                sep = " \n"
            )
        }
    }
    ## Constrain correlations
    if (XVar == TRUE & YVar == TRUE) {
        if (trait == TRUE) {
            modelConstraints <- paste(
                modelConstraints,
                "NEW(cor_txty);",
                "cor_txty = cov_txty / (sqrt(tx) * sqrt(ty)); \n",
                sep = " \n"
            )
        } else {
            modelConstraint <- paste(
                modelConstraints,
                "cov_txty = 0; \n",
                sep = " \n"
            )
        }
        if (AR == TRUE) {
            modelConstraints <- paste(
                modelConstraints,
                "NEW(cor_arxary);",
                "cor_arxary = cov_xy/(sqrt(arvx) * sqrt(arvy)); \n",
                sep = " \n"
            )
        } else {
            modelConstraints <- paste(
                modelConstraints,
                "cov_xy = 0; \n",
                sep = " \n"
            )
        }
        if (state == TRUE) {
            modelConstraints <- paste(
                modelConstraints,
                "NEW(cor_sxsy);",
                "cor_sxsy = cov_s/(sqrt(sx) * sqrt(sy)); \n",
                sep = " \n"
            )
        } 
    }
    if (stationarity == TRUE | AR == FALSE) {
        ## X AR Variance
        if (XVar == TRUE) {
            modelConstraints <- paste(
                modelConstraints,
                "arvx2 = arvx - arvx*a*a - 2*a*d*cov_xy; \n",
                sep = " \n"
            )
            for (w in wavesList[-c(1,2)]) {
                modelConstraints <- paste0(
                    modelConstraints,
                    "arvx",
                    w,
                    " = arvx2; \n"
                )
            }
        }
        ## Y AR Variance
        if (YVar == TRUE) {
            modelConstraints <- paste(
                modelConstraints,
                "arvy2 = arvy - arvy*b*b - 2*b*c*cov_xy; \n",
                sep = " \n"
            )
            for (w in wavesList[-c(1,2)]) {
                modelConstraints <- paste0(
                    modelConstraints,
                    "arvy",
                    w,
                    " = arvy2; \n"
                )
            }
        } 
    }
    if (stationarity == TRUE | state == FALSE) {
        ## X state Variance
        if (XVar == TRUE) {
            modelConstraints <- paste0(modelConstraints, " \n")
            for (w in wavesList[-c(1)]) {
                modelConstraints <- paste0(
                    modelConstraints,
                    "sx",
                    w,
                    " = sx; \n"
                )
            }
        }
        ## Y AR Variance
        if (YVar == TRUE) {
            modelConstraints <- paste0(modelConstraints, " \n")
            for (w in wavesList[-c(1)]) {
                modelConstraints <- paste0(
                    modelConstraints,
                    "sy",
                    w,
                    " = sy; \n"
                )
            }
        } 
    }
    ## Covariances
    if (XVar == TRUE &
        YVar == TRUE &
        ((AR == TRUE &
          stationarity == TRUE) |
         AR == FALSE)) {
        modelConstraints <- paste(
            modelConstraints,
            "cov_xyr2 = (1-a*b-c*d)*cov_xy - a*c*arvx - b*d*arvy; \n",
            sep = " \n"
        )
        for (w in wavesList[-c(1, 2)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "cov_xyr",
                w,
                " = cov_xyr2; \n"
            )
        }
    }
    if (XVar == TRUE &
        YVar == TRUE &
        ((state == TRUE &
          stationarity == TRUE) |
         state == FALSE)) {
        modelConstraints <- paste0(modelConstraints, " \n")
        for (w in wavesList[-c(1)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "cov_s",
                w,
                " = cov_s; \n"
            )
        }
    }
    if (limits == TRUE) {
        modelConstraints <-
            paste(
                modelConstraints,
                buildLimits(
                    XVar = XVar,
                    YVar = YVar,
                    trait = trait,
                    AR = AR,
                    state = state,
                    stationarity = stationarity,
                    constrainCor = stationarity
                ),
                sep = " \n"
            )
    }
    if (XVar == TRUE) {
        modelConstraints <- paste0(
            modelConstraints,
            "\n",
            "NEW(trait.p.x); \n",
            "NEW(ar.p.x); \n",
            "NEW(st.p.x); \n",
            "NEW(stab.x); \n",
            "trait.p.x = tx/(tx + arvx + sx); \n",
            "ar.p.x = arvx/(tx + arvx + sx); \n",
            "st.p.x = sx/(tx + arvx + sx); \n",
            "stab.x = a; \n"
        )
    }
    if (YVar == TRUE) {
        modelConstraints <- paste0(
            modelConstraints,
            "\n",
            "NEW(trait.p.y); \n",
            "NEW(ar.p.y); \n",
            "NEW(st.p.y); \n",
            "NEW(stab.y); \n",
            "trait.p.y = ty/(ty + arvy + sy); \n",
            "ar.p.y = arvy/(ty + arvy + sy); \n",
            "st.p.y = sy/(ty + arvy + sy); \n",
            "stab.y = b; \n"
        )
    }
    return(modelConstraints)
}


buildLimits <- function(XVar = TRUE,
                        YVar = TRUE,
                        trait = TRUE,
                        AR = TRUE,
                        state = TRUE,
                        stationarity = TRUE,
                        constrainCor = stationarity) {
    limits <- " \n"
    if (trait == TRUE & XVar == TRUE) {
        limits <- paste(
            limits,
            "tx > 0;",
            sep = " \n"
        )
    }
    if (trait == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "ty > 0;",
            sep = " \n"
        )
    }
    if (AR == TRUE & XVar == TRUE) {
        limits <- paste(
            limits,
            "arvx > 0;",
            sep = " \n"
        )
    }
    if (AR == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "arvy > 0;",
            sep = " \n"
        )
    }
    if (state == TRUE & XVar == TRUE) {
        limits <- paste(
            limits,
            "sx > 0;",
            sep = " \n"
        )
    }
    if (state == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "sy > 0;",
            sep = " \n"
        )
    }
    if (trait == TRUE & XVar == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "cor_txty < 1;",
            "cor_txty > -1;",
            sep = " \n"
        )
    }
    if (AR == TRUE & XVar == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "cor_arxary < 1;",
            "cor_arxary > -1;",
            sep = " \n"
        )
    }
    if (state == TRUE & XVar == TRUE & YVar == TRUE & stationarity == TRUE) {
        limits <- paste(
            limits,
            "cor_sxsy < 1;",
            "cor_sxsy > -1;",
            sep = " \n"
        )
    }
    return(limits)
}


#' Runs Mplus Models
#'
#' `run_starts_mplus()` produces code for variations of the general Stable Trait,
#' Autoregressive Trait, State Model and then runs the model using
#' MplusAutomation. Univariate and bivariate versions can be specified. 
#' Reduced versions of the model can be specified by setting variance components
#' for trait, AR, or state to be `FALSE`. 
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param XVar Logical. Include X variable.
#' @param YVar logical. Include Y variable.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param trait Logical value indicating whether to include a stable trait.
#'   Defaults to TRUE.
#' @param AR logical value indicating whether to include autoregressive trait.
#'   Defaults to TRUE.
#' @param state logical value indicating whether to include state. Defaults to
#'   TRUE.
#' @param crossLag logical value indicating whether to include cross-lagged
#'   paths. Defaults to TRUE.
#' @param stateCor logical value indicating whether to include state
#'   correlations. Defaults to TRUE.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCor logical value indicating whether to constrain state
#'   correlations to be the same.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `.#buildMplus.R.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_starts_mplus <- function(data,
                             waves,
                             XVar = TRUE,
                             YVar = TRUE,
                             xWaves = NULL,
                             yWaves = NULL,
                             xIndicators = 1,
                             yIndicators = 1,
                             trait = TRUE,
                             AR = TRUE,
                             state = TRUE,
                             crossLag = TRUE,
                             stateCor = TRUE,
                             stationarity = TRUE,
                             constrainCor = TRUE,
                             limits = TRUE,
                             dir="mplus",
                             title="test",
                             analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                             output = "stdyx; \n  cinterval; \n") {
    ## Set file name for output
    outFile <- paste(dir, title, sep="/")

    ## Build Model and Constraints
    modelStatement <- buildMplus(waves = waves,
                                 xWaves = xWaves,
                                 yWaves = yWaves,
                                 xIndicators = xIndicators,
                                 yIndicators = yIndicators,
                                 trait = trait,
                                 XVar = XVar,
                                 YVar = YVar,
                                 crossLag = crossLag,
                                 state = state,
                                 AR = AR,
                                 stateCor = stateCor,
                                 stationarity = stationarity)
    constraints <- buildConstraints(waves = waves,
                                    XVar = XVar,
                                    YVar = YVar,
                                    crossLag = crossLag,
                                    trait = trait,
                                    AR = AR,
                                    state = state,
                                    stationarity = stationarity,
                                    constrainCor = constrainCor,
                                    limits = limits)
    ## Clarify which waves exist
    if (is.null(xWaves)) xWaves <- 1:waves
    if (is.null(yWaves)) yWaves <- 1:waves
    ## Set labels for indicators
    indLabelsX <- letters[1:xIndicators]
    indLabelsY <- letters[1:yIndicators]
    if (XVar == TRUE) {
        if (xIndicators == 1) {
            xVariables <- paste0("x", xWaves)
        } else {
            xVariables <- paste0("x",
                                 rep(xWaves, each=xIndicators),
                                 indLabelsX)
        }
    } else {
        xVariables <- NULL
    }
    
    if (YVar == TRUE) {
        if (yIndicators == 1) {
            yVariables <- paste0("y", yWaves)
        } else {
            yVariables <- paste0("y",
                                 rep(yWaves, each=yIndicators),
                                 indLabelsY)
        }
    } else {
        yVariables <- NULL
    }
    variables <- c(xVariables, yVariables)
        
    inp <- MplusAutomation::mplusObject(TITLE=title,
                       rdata=data,
                       usevariables = variables,
                       ANALYSIS = analysis,
                       MODEL=modelStatement,
                       MODELCONSTRAINT = constraints,
                       OUTPUT=output
                       )
    output <- MplusAutomation::mplusModeler(inp, modelout = paste0(outFile, ".inp"), run=1)
    modelOutput <- MplusAutomation::readModels(target = paste0(dir, "/", title, ".out"))
    print("Errors:")
    print(output$Errors)
    print("Warnings:")
    print(output$Warnings)
    ##params <- modelOutput$parameters$unstandardized
    ##print(params)
    return(output)
}

#' Runs RI-CLPM in Mplus
#'
#' `run_riclpm_mplus()` Wrapper function for `run_starts_mplus()`. Produces and
#' runs Mplus code for the Random-Intercept Cross-Lagged Panel Model, which is a
#' STARTS model where the state variance is set to 0.
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCor logical value indicating whether to constrain state
#'   correlations to be the same.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `mplus`.
#' @param title Character vector. Defaults to `riclpm`.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_riclpm_mplus <- function(data,
                             waves,
                             xWaves = NULL,
                             yWaves = NULL,
                             xIndicators = 1,
                             yIndicators = 1,
                             stationarity = TRUE,
                             constrainCor = TRUE,
                             limits = TRUE,
                             dir = "mplus",
                             title = "riclpm",
                             analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                             output = "stdyx; \n  cinterval; \n"
                             ) {
    run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        state = FALSE,
        stationarity = stationarity,
        constrainCor = constrainCor,
        limits = limits,
        dir = dir,
        title = title,
        analysis = analysis,
        output = output
    )
}


#' Runs ARTS Model in Mplus
#'
#' `run_arts_mplus()` Wrapper function for `run_starts_mplus()`. Produces and
#' runs Mplus code for the ARTS Model, which is a STARTS model where the stable
#' trait variance is set to 0. This is also known as a factor CLPM.
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCor logical value indicating whether to constrain state
#'   correlations to be the same.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `mplus`.
#' @param title Character vector. Defaults to `arts`.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_arts_mplus <- function(data,
                           waves,
                           xWaves = NULL,
                           yWaves = NULL,
                           xIndicators = 1,
                           yIndicators = 1,
                           stationarity = TRUE,
                           constrainCor = TRUE,
                           limits = TRUE,
                           dir = "mplus",
                           title = "arts",
                           analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                           output = "stdyx; \n  cinterval; \n") {
    run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        stationarity = stationarity,
        constrainCor = constrainCor,
        trait = FALSE,
        state = TRUE,
        dir = dir,
        title = title,
        analysis = analysis,
        output = output
    )
}

#' Runs STS Model in Mplus
#'
#' `run_sts_mplus()` Wrapper function for `run_starts_mplus()`. Produces and
#' runs Mplus code for the STS Model, which is a STARTS model where the
#' autoregressive trait variance is set to 0. This is also known as a factor
#' CLPM.
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCor logical value indicating whether to constrain state
#'   correlations to be the same.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `mplus`.
#' @param title Character vector. Defaults to `sts`.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_sts_mplus <- function(data,
                          waves,
                          xWaves = NULL,
                          yWaves = NULL,
                          xIndicators = 1,
                          yIndicators = 1,
                          stationarity = TRUE,
                          constrainCor = TRUE,
                          limits = TRUE,
                          dir = "mplus",
                          title = "sts",
                          analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                          output = "stdyx; \n  cinterval; \n") {
    run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        trait = TRUE,
        state = TRUE,
        AR = FALSE,
        stationarity = stationarity,
        constrainCor = constrainCor,
        dir = dir,
        title = title,
        output = output,
        analysis = analysis
    )
}


#' Runs CLPM Model in Mplus
#'
#' `run_clpm_mplus()` Wrapper function for `run_starts_mplus()`. Produces and
#' runs Mplus code for the CLPM Model, which is a STARTS model where stable
#' trait and state variance is set to 0. 
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCor logical value indicating whether to constrain state
#'   correlations to be the same.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `mplus`.
#' @param title Character vector. Defaults to `sts`.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_clpm_mplus <- function(data,
                           waves,
                           xWaves = NULL,
                           yWaves = NULL,
                           xIndicators = 1,
                           yIndicators = 1,
                           stationarity = TRUE,
                           constrainCor = TRUE,
                           limits = TRUE,
                           dir = "mplus",
                           title = "clpm",
                           analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                           output = "stdyx; \n  cinterval; \n") {
    run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        trait = FALSE,
        state = FALSE,
        AR = TRUE,
        stationarity = stationarity,
        constrainCor = constrainCor,
        limits = limits,
        dir = dir,
        title = title,
        analysis = analysis,
        output = output
    )
}

#' Runs Univariate Starts Model for X in Mplus
#'
#' `run_sts_mplus()` Wrapper function for `run_starts_mplus()`. Produces and
#' runs Mplus code for the STS Model, which is a STARTS model where the
#' autoregressive trait variance is set to 0. This is also known as a factor
#' CLPM.
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `mplus`.
#' @param title Character vector. Defaults to `sts`.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_startsx_mplus <- function(data,
                              waves,
                              xWaves = NULL,
                              xIndicators = 1,
                              stationarity = TRUE,
                              limits = TRUE,
                              dir = "mplus",
                              title = "startsx",
                              analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                              output = "stdyx; \n  cinterval; \n") {
    run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        xIndicators = xIndicators,
        trait = TRUE,
        state = TRUE,
        YVar = FALSE,
        AR = TRUE,
        stationarity = stationarity,
        limits = limits,
        dir = dir,
        title = title,
        analysis = analysis,
        output = output
    )
}

#' Runs Univariate Starts Model for X in Mplus
#'
#' `run_sts_mplus()` Wrapper function for `run_starts_mplus()`. Produces and
#' runs Mplus code for the STS Model, which is a STARTS model where the
#' autoregressive trait variance is set to 0. This is also known as a factor
#' CLPM.
#'
#' @param data Dataframe that contains variables for analysis.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `mplus`.
#' @param title Character vector. Defaults to `sts`.
#' @param analysis Character vector listing the analysis statement for Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `MODEL=NOCOVARIANCES;\nCOVERAGE=.001;` .
#' @param output Character vector listing options for the output line in Mplus.
#'   Separate options with semicolon and `\n`. Defaults to
#'   `stdyx; \n  cinterval; \n`.
#' @returns Returns output from MplusAutomation::MplusModeler command, including
#'   estimates and fit statistics.
#' @export
run_startsy_mplus <- function(data,
                              waves,
                              yWaves = NULL,
                              yIndicators = 1,
                              stationarity = TRUE,
                              limits = TRUE,
                              dir = "mplus",
                              title = "startsy",
                              analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                              output = "stdyx; \n  cinterval; \n") {
    run_starts_mplus(
        data = data,
        waves = waves,
        yWaves = yWaves,
        yIndicators = yIndicators,
        trait = TRUE,
        state = TRUE,
        XVar = FALSE,
        AR = TRUE,
        stationarity = stationarity,
        limits = limits,
        dir = dir,
        title = title,
        analysis = analysis,
        output = output
    )
}

