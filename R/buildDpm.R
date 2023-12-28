################################################################################
## Function to create and run mplus DPM models
## Can include phantom variables for missing waves
## 
## "waves" = total number of possible waves (e.g., 20 for HILDA, 37 for SOEP)
## "actual waves" = waves with data (e.g., c(1,2,3,5,7,9))
################################################################################

################################################################################
## Build trait part of input file
################################################################################
buildTraitX_d <- function(waves) {
    wavesList <- 2:waves
    title <- "!!! Stable Trait X;\n"
    ## First do "BY" statements
    traitModel <- paste0(title,"traitX by arx", wavesList[1],"@1;\n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(traitModel,
                         paste0("traitX by arx", w, "@1;\n")
                         )
    }
    ## Then set variance
    traitModel <- paste0(traitModel, "\n!!! Trait Variance for X \ntraitX (tx);\n")
    return(traitModel)
}

buildTraitY_d <- function(waves) {
    wavesList <- 2:waves
    title <- "!!! Stable Trait Y;\n"
    ## First do "BY" statements
    traitModel <- paste0(title,"traitY by ary", wavesList[1],"@1;\n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(traitModel,
                         paste0("traitY by ary", w, "@1;\n")
                         )
    }
    ## Then set variance
    traitModel <- paste0(traitModel, "\n!!! Trait Variance for Y \ntraitY (ty);\n")
    return(traitModel)
}


################################################################################
## Build autoregressive part of input file
################################################################################

buildARX_d <- function(waves) {
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


buildARY_d <- function(waves) {
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
## Build observed part of model
################################################################################

buildObservedX_d <- function(xWaves, xIndicators=1) {
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

buildObservedY_d <- function(yWaves, yIndicators=1) {
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

buildStateX_d <- function(waves) {
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

buildStateY_d <- function(waves) {
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
buildLatentVarX_d <- function(waves) {
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
buildLatentVarY_d <- function(waves) {
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

buildClYFromX_d <- function(wavesList) {
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


buildClXFromY_d <- function(wavesList) {
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

buildTraitCors_d <- function() {
    traitCors <- "\n!!! Stable Trait Correlations \ntraitX with traitY(cov_txty); \n"
    return(traitCors)
}

buildArCors_d <- function(wavesList) {
    arCors <- "!!! AR Correlations \narx1 with ary1(cov_xy); \n"
    for (w in wavesList[-1]) {
        arCors <- paste0(
            arCors,
            paste0(
                "arx",
                w,
                " with ary",
                w,
                "; \n"
            )
        )
    }
    return(arCors)
}

buildTraitArCors_d <- function(wavesList, XVar=TRUE, YVar=TRUE) {
    cors <- "!!! Correlations between Stable Traits and AR; \n"
    if (XVar == TRUE) {
        cors <- paste0(
            cors,
            "traitX with arx1; \n"
        )
    }
    if (YVar == TRUE) {
        cors <- paste0(
            cors,
            "traitY with ary1; \n"
        )
    }
    if (XVar == TRUE & YVar == TRUE) {
        cors <- paste0(
            cors,
            "traitX with ary1; \n",
            "traitY with arx1; \n"
        )
    }
}
    
buildStateCors_d <- function(wavesList) {
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


buildResidCorsX_d <- function(xWaves, xIndicators, constrainCors) {
    residCors <- "!!! Residual Correlations \n"
    indLabels <- letters[1:xIndicators]
    if (constrainCors == TRUE) {
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
    } else {
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
                        "; \n"
                    )
                }
            }
        }
    }
    return(residCors)
}

buildResidCorsY_d <- function(yWaves, yIndicators, constrainCors = TRUE) {
    residCors <- "!!! Residual Correlations \n"
    indLabels <- letters[1:yIndicators]
    if (constrainCors == TRUE) {
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
                        "; \n"
                    )
                }
            }
        }
    } else {
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
                        "); \n"
                    )
                }
            }
        }
    }
    return(residCors)
}



################################################################################
## Build final MODEL statement
################################################################################
#' Builds Mplus Code for Dynamic Panel Model
#'
#' `buildDpm()` produces code for variations of the Dynamic Panel Model.
#' Univariate and bivariate versions can be specified. Cross-lagged paths are
#' included by default, though they can be removed. Correlations between state
#' and autoregressive components can also be included.
#'
#' @param waves Numeric value indicating total number of waves.
#' @param XVar Logical. Include X variable.
#' @param YVar logical. Include Y variable.
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param crossLag logical value indicating whether to include cross-lagged
#'   paths. Defaults to TRUE.
#' @param arCor logical value indicating whether to include autoregressive trait
#'   correlations. Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
buildDpm <- function(waves,
                     XVar = TRUE,
                     YVar = TRUE,
                     xIndicators = 1,
                     yIndicators = 1,
                     crossLag = TRUE,
                     arCor = TRUE,
                     constrainCors = TRUE) {
    xWaves <- 1:waves
    yWaves <- 1:waves
    wavesList <- 1:waves
    if (XVar == FALSE | YVar == FALSE) {
        crossLag <- FALSE
        arCor <- FALSE
    }
    ## Set state variance to false; consider implementing in future version.
    state <- FALSE
    ## Start model statement
    modelStatement <- "!!! Automated Code for Mplus \n"

    ## Add latent variables and indicators
    if (XVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildObservedX_d(xWaves, xIndicators),
            buildLatentVarX_d(waves),
            buildARX_d(waves),
            buildTraitX_d(waves),
            buildStateX_d(waves),
            sep = " \n"
        )
    }
    if (YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildObservedY_d(yWaves, yIndicators),
            buildLatentVarY_d(waves),
            buildARY_d(waves),
            buildTraitY_d(waves),
            buildStateY_d(waves),
            sep = " \n"
        )
    }

    ## Add Cross Lags
    if (crossLag == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildClYFromX_d(wavesList),
            buildClXFromY_d(wavesList),
            sep = " \n"
        )
    }

    ## Add Trait Cors
    if (XVar == TRUE & YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildTraitCors_d(),
            sep = " \n"
        )
    }

    ## Add AR Cors
    if (XVar == TRUE & YVar == TRUE & arCor == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildArCors_d(wavesList),
            sep = " \n"
        )
    }

    ## Add Trait AR Cors
    modelStatement <- paste(
        modelStatement,
        buildTraitArCors_d(wavesList, XVar, YVar)
    )

    ## Add State Cors; later set to 0 if needed
    if (XVar == TRUE & YVar == TRUE) {
        modelStatement <- paste(
            modelStatement,
            buildStateCors_d(wavesList),
            sep = " \n"
        )
    }

    ## Add Resid Cors
    if (XVar == TRUE & xIndicators > 1) {
        modelStatement <- paste(
            modelStatement,
            buildResidCorsX_d(xWaves, xIndicators, constrainCors),
            sep = " \n"
        )
    }

    if (YVar == TRUE & yIndicators > 1) {
        modelStatement <- paste(
            modelStatement,
            buildResidCorsY_d(yWaves, yIndicators, constrainCors),
            sep = " \n"
        )
    }

    ## Return model
    return(modelStatement)
}


## Function to specify constraints

#' Build Constraints for DPM
#'
#' Builds constraint statement for Mplus input file based on specified options.
#' 
#' @param waves Numeric value indicating total number of waves.
#' @param XVar Logical. Include X variable.
#' @param YVar logical. Include Y variable.
#' @param crossLag logical value indicating whether to include cross-lagged
#'   paths. Defaults to TRUE.
#' @param arCor Logical value indicating whether to include within-wave
#'   correlations for AR traits.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to valid values. Defaults to TRUE.
#' @returns Character vector representing the Mplus code for the constraints
#'   statement.
#' @export
buildConstraints_d <- function(waves,
                             XVar = TRUE,
                             YVar = TRUE,
                             crossLag = TRUE,
                             arCor = TRUE,
                             limits = TRUE) {
    wavesList <- 1:waves
    if (XVar == FALSE | YVar == FALSE) {
        crossLag <- FALSE
    }
    ## Consider implementing state variance; disabled for now.
    state <- FALSE
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
    if (XVar == TRUE) {
        for (w in wavesList[-c(1, waves)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "a",
                w,
                " = a; \n"
            )
        }
    }
    if (YVar == TRUE) {
        for (w in wavesList[-c(1, waves)]) {
            modelConstraints <- paste0(
                modelConstraints,
                "b",
                w,
                " = b; \n"
            )
        }
    }
    if (crossLag == TRUE) {
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
    modelConstraints <- paste(
        modelConstraints,
        "NEW(cor_arxary);",
        "cor_arxary = cov_xy/(sqrt(arvx) * sqrt(arvy)); \n",
        sep = " \n"
    )
    if (state == TRUE) {
            modelConstraints <- paste(
                modelConstraints,
                "NEW(cor_sxsy);",
                "cor_sxsy = cov_s/(sqrt(sx) * sqrt(sy)); \n",
                sep = " \n"
            )
    }
    if (state == FALSE) {
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
        ## Y state Variance
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
    if (XVar == TRUE &
        YVar == TRUE) {
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
                buildLimits_d(
                    XVar = XVar,
                    YVar = YVar
                ),
                sep = " \n"
            )
    }
    if (XVar == TRUE) {
        modelConstraints <- paste0(
            modelConstraints,
            "\n",
            "NEW(trait.x); \n",
            "NEW(ar.p.x); \n",
            "NEW(st.p.x); \n",
            "NEW(stab.x); \n",
            "trait.x = tx/(tx + arvx + sx); \n",
            "ar.p.x = arvx/(tx + arvx + sx); \n",
            "st.p.x = sx/(tx + arvx + sx); \n",
            "stab.x = a; \n"
        )
    }
    if (YVar == TRUE) {
        modelConstraints <- paste0(
            modelConstraints,
            "\n",
            "NEW(trait.y); \n",
            "NEW(ar.p.y); \n",
            "NEW(st.p.y); \n",
            "NEW(stab.y); \n",
            "trait.y = ty/(ty + arvy + sy); \n",
            "ar.p.y = arvy/(ty + arvy + sy); \n",
            "st.p.y = sy/(ty + arvy + sy); \n",
            "stab.y = b; \n"
        )
    }
    return(modelConstraints)
}


buildLimits_d <- function(XVar = TRUE,
                          YVar = TRUE) {
    state <- FALSE
    limits <- " \n"
    if (XVar == TRUE) {
        limits <- paste(
            limits,
            "tx > 0;",
            sep = " \n"
        )
    }
    if (YVar == TRUE) {
        limits <- paste(
            limits,
            "ty > 0;",
            sep = " \n"
        )
    }
    if (XVar == TRUE) {
        limits <- paste(
            limits,
            "arvx > 0;",
            sep = " \n"
        )
    }
    if (YVar == TRUE) {
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
    if (XVar == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "cor_txty < 1;",
            "cor_txty > -1;",
            sep = " \n"
        )
    }
    if (XVar == TRUE & YVar == TRUE) {
        limits <- paste(
            limits,
            "cor_arxary < 1;",
            "cor_arxary > -1;",
            sep = " \n"
        )
    }
    if (state == TRUE & XVar == TRUE & YVar == TRUE) {
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
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1.
#' @param crossLag logical value indicating whether to include cross-lagged
#'   paths. Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @param dir Character vector listing directory for mplus files. Defaults to
#'   `.#buildMplus.R.
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
run_dpm_mplus <- function(data,
                          waves,
                          XVar = TRUE,
                          YVar = TRUE,
                          xIndicators = 1,
                          yIndicators = 1,
                          crossLag = TRUE,
                          constrainCors = TRUE,
                          limits = TRUE,
                          dir = "mplus",
                          title = "dpm",
                          analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;",
                          output = "stdyx; \n  cinterval; \n") {
    ## Set file name for output
    outFile <- paste(dir, title, sep = "/")

    ## Build Model and Constraints
    modelStatement <- buildDpm(
        waves = waves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        XVar = XVar,
        YVar = YVar,
        crossLag = crossLag,
        constrainCors = constrainCors
    )
    constraints <- buildConstraints_d(
        waves = waves,
        XVar = XVar,
        YVar = YVar,
        crossLag = crossLag,
        limits = limits
    )
    ## Clarify which waves exist
    xWaves <- 1:waves
    yWaves <- 1:waves
    ## Set labels for indicators
    indLabelsX <- letters[1:xIndicators]
    indLabelsY <- letters[1:yIndicators]
    if (XVar == TRUE) {
        if (xIndicators == 1) {
            xVariables <- paste0("x", xWaves)
        } else {
            xVariables <- paste0(
                "x",
                rep(xWaves, each = xIndicators),
                indLabelsX
            )
        }
    } else {
        xVariables <- NULL
    }

    if (YVar == TRUE) {
        if (yIndicators == 1) {
            yVariables <- paste0("y", yWaves)
        } else {
            yVariables <- paste0(
                "y",
                rep(yWaves, each = yIndicators),
                indLabelsY
            )
        }
    } else {
        yVariables <- NULL
    }
    variables <- c(xVariables, yVariables)

    inp <- MplusAutomation::mplusObject(
        TITLE = title,
        rdata = data,
        usevariables = variables,
        ANALYSIS = analysis,
        MODEL = modelStatement,
        MODELCONSTRAINT = constraints,
        OUTPUT = output
    )

    output <- MplusAutomation::mplusModeler(inp, modelout = paste0(outFile, ".inp"), run=1)
    outputList <- list(mplusOutput = output,
                       xTrait = output$results$parameters$unstandardized[
                                                              which(output$results$parameters$unstandardized$param == "TRAIT.X"),
                                                              "est"
                                                          ],
                       xAr = output$results$parameters$unstandardized[
                                                           which(output$results$parameters$unstandardized$param == "AR.P.X"),
                                                           "est"
                                                       ],
                       xStability = output$results$parameters$unstandardized[
                                                                  which(output$results$parameters$unstandardized$param == "STAB.X"),
                                                                  "est"
                                                              ],
                       yTrait = output$results$parameters$unstandardized[
                                                              which(output$results$parameters$unstandardized$param == "TRAIT.Y"),
                                                              "est"
                                                          ],
                       yAr = output$results$parameters$unstandardized[
                                                           which(output$results$parameters$unstandardized$param == "AR.P.Y"),
                                                           "est"
                                                       ],
                       yStability = output$results$parameters$unstandardized[
                                                                  which(output$results$parameters$unstandardized$param == "STAB.Y"),
                                                                  "est"
                                                              ],
                       traitCor = output$results$parameters$unstandardized[
                                                                which(output$results$parameters$unstandardized$param == "COR_TXTY"),
                                                                "est"
                                                            ],
                       arCor = output$results$parameters$unstandardized[
                                                             which(output$results$parameters$unstandardized$param == "COR_ARXA"),
                                                             "est"
                                                         ],
                       tXarX = output$results$parameters$stdyx.standardized[
                                                             which(output$results$parameters$stdyx.standardized$paramHeader == "TRAITX.WITH" &
                                                                   output$results$parameters$stdyx.standardized$param == "ARX1"),
                                                             "est"
                                                         ],
                       tXarY = output$results$parameters$stdyx.standardized[
                                                             which(output$results$parameters$stdyx.standardized$paramHeader == "TRAITX.WITH" &
                                                                   output$results$parameters$stdyx.standardized$param == "ARY1"),
                                                             "est"
                                                         ],
                       tYarY = output$results$parameters$stdyx.standardized[
                                                             which(output$results$parameters$stdyx.standardized$paramHeader == "TRAITY.WITH" & output$results$parameters$stdyx.standardized$param == "ARY1"),
                                                             "est"
                                                         ],
                       tYarX = output$results$parameters$stdyx.standardized[
                                                             which(output$results$parameters$stdyx.standardized$paramHeader == "TRAITY.WITH" &
                                                                   output$results$parameters$stdyx.standardized$param == "ARX1"),
                                                             "est"
                                                         ],
                       yOnX = output$results$parameters$unstandardized[
                                                                which(output$results$parameters$unstandardized$paramHeader == "ARY2.ON" &
output$results$parameters$unstandardized$param == "ARX1"), 
                                                                "est"
                                                            ],
                       xOnY = output$results$parameters$unstandardized[
                                                                which(output$results$parameters$unstandardized$paramHeader == "ARX2.ON" &
output$results$parameters$unstandardized$param == "ARY1"), 
                                                                "est"
                                                            ]
                       )
    class(outputList) <- "pcmdObject"
    print("Errors:")
    print(output$Errors)
    print("Warnings:")
    print(output$Warnings)
    summary(outputList)
    return(outputList)
}


#' Summarizes results from `run_dpm_mplus()`
#'
#' @param object Results from `run_dpm_mplus()`.
#' @param ... Additional arguments to `summary()`.
#'
#' @export
summary.pcmdObject <- function(object, ...) {
    cat("Model Summary: \n")
    summary(object$mplusOutput)
    cat("\n")
    cat("Variance Decomposition (First Wave): \n")
    cat("\n")
    cat("X Variable: \n")
    cat("Trait: ", object$xTrait,
        ", Autoregressive: ",object$xAr, "\n")
    cat("Y Variable: \n")
    cat("Trait: ", object$yTrait,
        ", Autoregressive: ", object$yAr, "\n")
    cat("\n")
    cat("Stability: \n")
    cat("X: ", object$xStab, "\n")
    cat("Y: ", object$yStab, "\n")
    cat("\n")
    cat("Cross-Lag Paths: \n")
    cat("Y predicted from X: ", object$yOnX ,"\n")
    cat("X predicted from Y: ", object$xOnY,"\n")
    cat("\n")
    cat("Correlations: \n")
    cat("Stable Trait: ", object$traitCor, "\n")
    cat("Autoregressive Trait: ", object$arCor, "\n")
    cat("X Trait with X AR: ", object$tXarX, "\n")
    cat("X Trait with Y AR: ", object$tXarY, "\n")
    cat("Y Trait with Y AR: ", object$tYarY, "\n")
    cat("Y Trait with X AR: ", object$tYarX, "\n")
    }


    
