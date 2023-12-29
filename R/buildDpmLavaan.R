################################################################################
## Function to create and run lavaan DPM models
## 
## "waves" = total number of possible waves (e.g., 20 for HILDA, 37 for SOEP)
## "actual waves" = waves with data (e.g., c(1,2,3,5,7,9))
################################################################################

################################################################################
## Model building starts here
################################################################################

## Build trait part of input file
buildTraitX_ld <- function(wavesList) {
    
    traitModel <- "## Stable Trait/Random Intercepts for X \n"
    ## First do indicator statements
    for (w in wavesList[-1]) {
        traitModel <- paste0(
            traitModel,
            paste0("ri_x =~ 1*arx", w, " \n")
        )
    }
    ## Then set variance
    traitModel <- paste0(
        traitModel,
        "## Trait Variance \nri_x ~~ x_t_var*ri_x \n" ## Constrained to be >0
    )
    return(traitModel)
}


buildTraitY_ld <- function(wavesList) {
    traitModel <- "## Stable Trait/Random Intercepts for Y \n"
    ## First do indicator statements
    for (w in wavesList[-1]) {
        traitModel <- paste0(
            traitModel,
            paste0("ri_y =~ 1*ary", w, " \n")
        )
    }
    ## Then set variance
    traitModel <- paste0(
        traitModel,
        "## Trait Variance \nri_y ~~ y_t_var*ri_y \n" ## Constrainted to be >0
    )
    return(traitModel)
}


################################################################################
## Build latent occasion/observed part of model
################################################################################

buildObservedX_ld <- function(xWaves, xIndicators) {
    title <- "### Latent Occasion Variables for X \n"
    ## First do indicator statements
    observedModel <- title
    if (xIndicators == 1) {
        for (w in xWaves) {
            observedModel <- paste0(
                observedModel,
                "lx",
                w,
                " =~ 1*x",
                w,
                "\n"
            )
        }
        ## Then constrain variances to 0
        observedModel <- paste0(
            observedModel,
            "## Residual Variance Constrained to 0 \n"
        )
        for (w in xWaves) {
            observedModel <- paste0(
                observedModel,
                "x",
                w,
                " ~~ ",
                0,
                "*x",
                w,
                " \n"
            )
        }
    } else {
        indLabels <- letters[1:xIndicators]
        for (w in xWaves) {
            for (i in indLabels) {
                if (i == "a") {
                    observedModel <- (
                        paste0(
                            observedModel,
                            "lx",
                            w,
                            " =~ 1*x",
                            w,
                            i,
                            "\n"
                        )
                    )
                } else {
                    observedModel <- (
                        paste0(
                            observedModel,
                            "lx",
                            w,
                            " =~ lx",
                            i,
                            "*x",
                            w,
                            i,
                            "\n"
                        )
                    )
                }
            }
        }
        observedModel <- paste0(
            observedModel,
            "## Residual Variance for Indicators \n"
        )
        for (w in xWaves) {
            for (i in indLabels) {
                observedModel <- paste0(
                    observedModel,
                    "x", w, i,
                    " ~~ ",
                    "xr", i,
                    "*x", w, i,
                    " \n"
                )
            }
        }
    }
    return(observedModel)
}

buildObservedY_ld <- function(yWaves, yIndicators) {
    title <- "### Latent Occasion Variables for Y \n"
    ## First do indicator statements
    observedModel <- title
    if (yIndicators == 1) {
        for (w in yWaves) {
            observedModel <- paste0(
                observedModel,
                "ly",
                w,
                " =~ 1*y",
                w,
                "\n"
            )
        }
        ## Then constrain variances to 0
        observedModel <- paste0(
            observedModel,
            "\n## Residual Variance Constrained to 0 \n"
        )
        for (w in yWaves) {
            observedModel <- paste0(
                observedModel,
                "y",
                w,
                " ~~ ",
                0,
                "*y",
                w,
                " \n"
            )
        }
    } else {
        indLabels <- letters[1:yIndicators]
        for (w in yWaves) {
            for (i in indLabels) {
                if (i == "a") {
                    observedModel <- (
                        paste0(
                            observedModel,
                            "ly",
                            w,
                            " =~ 1*y",
                            w,
                            i,
                            "\n"
                        )
                    )
                } else {
                    observedModel <- (
                        paste0(
                            observedModel,
                            "ly",
                            w,
                            " =~ ly",
                            i,
                            "*y",
                            w,
                            i,
                            "\n"
                        )
                    )
                }
            }
        }
        observedModel <- paste0(
            observedModel,
            "\n## Residual Variance Constrained to 0 \n"
        )
        for (w in yWaves) {
            for (i in indLabels) {
                observedModel <- paste0(
                    observedModel,
                    "y", w, i,
                    " ~~ ",
                    "yr", i,
                    "*y", w, i,
                    " \n"
                )
            }
        }
    }
    return(observedModel)
}



################################################################################
## Build State Variance
################################################################################

buildStateX_ld <- function(wavesList, stationarity = TRUE) {
    title <- "\n### X States\n"
    stateModel <- title
    ## First do indicator statements
    for (w in wavesList) {
        stateModel <- paste0(
            stateModel,
            paste0(
                "sx",
                w,
                " =~ 1*lx",
                w,
                " \n"
            )
        )
    }
    ## Then add variances 
    stateModel <- paste0(
        stateModel,
        "\n## State Variance\n"
    )
    if (stationarity == TRUE) {
        for (w in wavesList) {
            stateModel <- paste0(
                stateModel,
                paste0(
                    "sx",
                    w,
                    " ~~ sx*sx",
                    w,
                    " \n"
                )
            )
        }
    } else {
        for (w in wavesList) {
            stateModel <- paste0(
                stateModel,
                paste0(
                    "sx",
                    w,
                    " ~~ sx",
                    w,
                    " \n"
                )
            )
        }
    }
    return(stateModel)
}




buildStateY_ld <- function(wavesList, stationarity = TRUE) {
    title <- "\n### Y States\n"
    stateModel <- title
    ## First do indicator statements
    for (w in wavesList) {
        stateModel <- paste0(
            stateModel,
            paste0(
                "sy",
                w,
                " =~ 1*ly",
                w,
                " \n"
            )
        )
    }
    ## Then add variances (constrained to be equal if stationarity = TRUE)
    stateModel <- paste0(
        stateModel,
        "\n## State Variance\n"
    )
    if (stationarity == TRUE) {
        for (w in wavesList) {
            stateModel <- paste0(
                stateModel,
                paste0(
                    "sy",
                    w,
                    " ~~ sy*sy",
                    w,
                    " \n"
                )
            )
        }
    } else {
        for (w in wavesList) {
            stateModel <- paste0(
                stateModel,
                paste0(
                    "sy",
                    w,
                    " ~~ sy",
                    w,
                    " \n"
                )
            )
        }
    }
    return(stateModel)
}

################################################################################
## latent occasion variance
################################################################################

buildLatentVarX_ld <- function(waves) {
    wavesList <- 1:waves
    title <- "\n## Constrain Latent Occasion Residuals for X to 0\n"
    latentVar <- paste0(title, "lx1 ~~ 0*lx1 \n")
    for (w in wavesList[-1]) {
        latentVar <- paste0(
            latentVar,
            paste0(
                "lx",
                w,
                " ~~ 0*lx",
                w,
                " \n"
            )
        )
    }
    return(latentVar)
}

buildLatentVarY_ld <- function(waves) {
    wavesList <- 1:waves
    title <- "\n## Constrain Latent Occasion Residuals for Y to 0\n"
    latentVar <- paste0(title, "ly1 ~~ 0*ly1 \n")
    for (w in wavesList[-1]) {
        latentVar <- paste0(
            latentVar,
            paste0(
                "ly",
                w,
                " ~~ 0*ly",
                w,
                " \n"
            )
        )
    }
    return(latentVar)
}



################################################################################
## Build autoregressive part of model
################################################################################

buildARX_ld <- function(waves) {
    wavesList <- 1:waves
    title <- "### Autoregressive Part for X\n"
    subtitle1 <- "## Indicator Statements\n"
    ## First do X indicator statements linking AR to latent occasion variable
    arXModel <- paste0(
        title,
        subtitle1,
        "arx",
        wavesList[1],
        " =~ 1*lx",
        wavesList[1],
        " \n"
    )
    for (w in wavesList[-1]) {
        arXModel <- paste0(
            arXModel,
            "arx",
            w,
            " =~ 1*lx",
            w,
            " \n"
        )
    }
    ## Then add X stability paths
    subtitle2 <- "\n## Regression Statements\n"
    arXModel <- paste0(
        arXModel,
        subtitle2,
        "arx2 ~ a*arx1 \n"
    )
    for (w in wavesList[-c(1:2)]) {
        arXModel <- paste0(
            arXModel,
            paste0(
                "arx",
                w,
                " ~ a*arx",
                (w - 1),
                " \n"
            )
        )
    }
    ## Finally add AR variances (waves after the first are set to be equal)
    arXModel <- paste0(
        arXModel,
        "\n## Autoregressive Component Variance \narx1 ~~ arvx*arx1 \n"
    )
    for (w in wavesList[-1]) {
        arXModel <- paste0(
            arXModel,
            paste0(
                "arx",
                w,
                " ~~ arx",
                w,
                " \n"
            )
        )
    }
    return(arXModel)
}

buildARY_ld <- function(waves) {
    wavesList <- 1:waves
    title <- "### Autoregressive Part for Y\n"
    subtitle1 <- "## Indicator Statements\n"
    ## First do Y indicator statements linking AR to latent occasion variable
    arYModel <- paste0(
        title,
        subtitle1,
        "ary",
        wavesList[1],
        " =~ 1*ly",
        wavesList[1],
        " \n"
    )
    for (w in wavesList[-1]) {
        arYModel <- paste0(
            arYModel,
            "ary",
            w,
            " =~ 1*ly",
            w,
            " \n"
        )
    }
    ## Then add Y stability paths
    subtitle2 <- "\n## Regression Statements\n"
    arYModel <- paste0(
        arYModel,
        subtitle2,
        "ary2 ~ b*ary1 \n"
    )
    for (w in wavesList[-c(1:2)]) {
        arYModel <- paste0(
            arYModel,
            paste0(
                "ary",
                w,
                " ~ b*ary",
                (w - 1),
                " \n"
            )
        )
    }
    ## Finally add AR variances (waves after the first are set to be equal)
    arYModel <- paste0(
        arYModel,
        "\n## Autoregressive Component Variance \nary1 ~~ arvy*ary1 \n"
    )
    for (w in wavesList[-1]) {
        arYModel <- paste0(
            arYModel,
            paste0(
                "ary",
                w,
                " ~~ ary",
                w,
                " \n"
            )
        )
    }
    return(arYModel)
}


################################################################################
## Build Correlations
################################################################################

buildTraitCors_ld <- function() {
    traitCors <- "## Stable Trait Correlations \n ri_x ~~ cov_t*ri_y \n"
    return(traitCors)
}


buildArCors_ld <- function(wavesList) {
    arCors <- "## AR Correlations \narx1 ~~ cor_xy*ary1 \n"
    for (w in wavesList[-1]) {
        arCors <- paste0(
            arCors,
            paste0("arx", w, " ~~ ary", w, " \n")
        )
    }
    return(arCors)
}

buildTraitArCors_ld <- function(wavesList, XVar=TRUE, YVar=TRUE) {
    cors <- "## Correlations between Stable Traits and AR \n"
    if (XVar == TRUE) {
        cors <- paste0(
            cors,
            "ri_x ~~ arx1 \n"
        )
    }
    if (YVar == TRUE) {
        cors <- paste0(
            cors,
            "ri_y ~~ ary1 \n"
        )
    }
    if (XVar == TRUE & YVar == TRUE) {
        cors <- paste0(
            cors,
            "ri_x ~~ ary1 \n",
            "ri_y ~~ arx1 \n"
            )
    }
    return(cors)
}


buildStateCors_ld <- function(wavesList) {
    sCors <- "## State Correlations \n"
    for (w in wavesList) {
        sCors <- paste0(
            sCors,
            paste0("sx", w, " ~~ sy", w, " \n")
        )
    }
    return(sCors)
}

buildResidCorsX_ld <- function(xWaves, xIndicators, constrainCors = TRUE) {
    residCors <- "## Residual Correlations \n"
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
                        " ~~ x",
                        i,
                        "lag",
                        k - xWaves[j],
                        " * x",
                        k,
                        i,
                        " \n"
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
                        " ~~ x",
                        i,
                        "lag",
                        k - xWaves[j],
                        "* x",
                        k,
                        i,
                        " \n"
                    )
                }
            }
        }
    }
    return(residCors)
}


buildResidCorsY_ld <- function(yWaves, yIndicators, constrainCors = TRUE) {
    residCors <- "## Residual Correlations \n"
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
                        " ~~ y",
                        i,
                        "lag",
                        k - yWaves[j],
                        " * y",
                        k,
                        i,
                        " \n"
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
                        " ~~ y",
                        i,
                        "lag",
                        k - yWaves[j],
                        "* y",
                        k,
                        i,
                        " \n"
                    )
                }
            }
        }
    }
    return(residCors)
}

                        


################################################################################
## Build Cross-Lagged Paths
################################################################################

buildClYFromX_ld <- function(wavesList) {
    title <- "### Cross-Lagged Path \n"
    subtitle <- "## Y predicted from X \n"
    clYonX <- paste0(title, subtitle)
    for (w in wavesList[-1]) {
        clYonX <- paste0(
            clYonX, "ary", w, " ~ c*arx", (w - 1), " \n"
        )
    }
    return(clYonX)
}


buildClXFromY_ld <- function(wavesList, stationarity = TRUE) {
    title <- "### Cross-Lagged Path \n"
    subtitle <- "## X predicted from Y \n"
    clXonY <- paste0(title, subtitle)
    for (w in wavesList[-1]) {
        clXonY <- paste0(
            clXonY, "arx", w, " ~ d*ary", (w - 1), " \n"
        )
    }
    return(clXonY)
}


## cons_cor_cl <- "cor_xyr == cor_xy -
## cor_xy * (a * b) -
##     cor_xy * (c * d) -
##     a * c * arvx -
## b * d * arvy \n"

    


################################################################################
## Build Final Model
################################################################################

#' Build Dynamic Panel Model Lavaan Code
#'
#' `buildLavaanDpm()` produces lavaan code for variations of the Dynamic Panel
#' Model. Univariate and bivariate versions can be specified. 
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
#' @param arCor Logical value indicating whether to include autoregresive trait
#'   correlations. Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
buildLavaanDpm <- function(waves,
                        XVar = TRUE,
                        YVar = TRUE,
                        xIndicators = 1,
                        yIndicators = 1,
                        crossLag = TRUE,
                        arCor = TRUE,
                        constrainCors = TRUE,
                        limits = TRUE) {
    xWaves <- 1:waves
    yWaves <- 1:waves
    wavesList <- 1:waves
    stationarity <- TRUE ## Not used, but included to reuse old code
    if (XVar == FALSE | YVar == FALSE) {
        crossLag <- FALSE
        arCor <- FALSE
    }
    state <- FALSE ## Set this to false to use old code
    modelCode <- "### Model Code Generated by buildModel.R"
    if (XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildObservedX_ld(xWaves, xIndicators),
            buildLatentVarX_ld(waves),
            sep=" \n")
    }
    if (YVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildObservedY_ld(yWaves, yIndicators),
            buildLatentVarY_ld(waves),
            sep = " \n"
        )
    }
    if (XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildTraitX_ld(wavesList)
        )
    }
    if (YVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildTraitY_ld(wavesList)
        )
    }
    if (state == TRUE & XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildStateX_ld(wavesList, stationarity),
            sep = " \n"
        )
    }
    if (state == TRUE & YVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildStateY_ld(wavesList, stationarity),
            sep = " \n"
        )
    }
    if (XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildARX_ld(waves),
            sep = " \n"
        )
    }
    if (YVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildARY_ld(waves),
            sep = " \n"
        )
    }
    if (YVar == TRUE  & XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildClYFromX_ld(wavesList),
            buildClXFromY_ld(wavesList),
            sep = " \n"
        )
    }
    if (YVar == TRUE & XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildTraitCors_ld(),
            sep = " \n"
        )
    }
    if (YVar == TRUE & XVar == TRUE) {
        modelCode <- paste(
            modelCode,
            buildArCors_ld(wavesList),
            sep = " \n"
        )
    }

    modelCode <- paste(
        modelCode,
        buildTraitArCors_ld(wavesList, XVar, YVar),
        sep = " \n"
    )
    
    if (XVar == TRUE &
        xIndicators > 1) {
        modelCode <- paste(
            modelCode,
            buildResidCorsX_ld(xWaves, xIndicators, constrainCors),
            sep = " \n"
        )
    }

    if (YVar == TRUE &
        yIndicators > 1) {
        modelCode <- paste(
            modelCode,
            buildResidCorsY_ld(yWaves, yIndicators, constrainCors),
            sep = " \n"
        )
    }

    
    ## ## Constraints
    ## cons_uni <- "arv2x == arvx - (a^2) * arvx \n"
    ## cons_y_uni <- "arv2y == arvy - (b^2) * arvy \n"
    ## cons_cl <- "arv2x == (1-a^2)*arvx - 2*a*d*cor_xy - d^2*arvy \n"
    ## cons_y_cl <- "arv2y == (1-b^2)*arvy - 2*b*c*cor_xy - c^2*arvx \n"
    ## cons_cor_noCl <- "cor_xyr == cor_xy - (a*cor_xy*b) \n"
    ## cons_cor_cl <- "cor_xyr == (1-a*b-c*d)*cor_xy - a*c*arvx - b*d*arvy \n"


    modelCode <- paste(
        modelCode,
        "## Constraints",
        sep = " \n"
    )

    ## Restrict variance to be non-zero
    if (limits == TRUE) {
        if (XVar == TRUE) {
            modelCode <- paste(
                modelCode,
                "x_t_var > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE) {
            modelCode <- paste(
                modelCode,
                "y_t_var > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE) {
            modelCode <- paste(
                modelCode,
                "arvx > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE) {
            modelCode <- paste(
                modelCode,
                "arvy > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE &
            state == TRUE) {
            modelCode <- paste(
                modelCode,
                "sx > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE &
            state == TRUE) {
            modelCode <- paste(
                modelCode,
                "sy > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE &
            YVar == TRUE &
            arCor == TRUE) {
            modelCode <- paste(
                modelCode,
                "cov_t < .99*(sqrt(x_t_var)*sqrt(y_t_var))",
                "cov_t > -.99*(sqrt(x_t_var)*sqrt(y_t_var))",
                sep = " \n"
            )
            modelCode <- paste(
                modelCode,
                "cor_xy < .99*(sqrt(arvx)*sqrt(arvy))",
                "cor_xy > -.99*(sqrt(arvx)*sqrt(arvy))",
                sep = " \n"
            )
        }
    }
    return(modelCode)
}

