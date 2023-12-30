################################################################################
## Function to create and run lavaan STARTS models
## 
## "waves" = total number of possible waves (e.g., 20 for HILDA, 37 for SOEP)
## "actual waves" = waves with data (e.g., c(1,2,3,5,7,9))
################################################################################

################################################################################
## Model building starts here
################################################################################


## Build trait part of input file
buildTraitX_l <- function(wavesList) {
    title <- "## Stable Trait/Random Intercepts for X \n"
    ## First do indicator statements
    traitModel <- paste0(title, "ri_x =~ 1*lx", wavesList[1], " \n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(
            traitModel,
            paste0("ri_x =~ 1*lx", w, " \n")
        )
    }
    ## Then set variance
    traitModel <- paste0(
        traitModel,
        "## Trait Variance \nri_x ~~ x_t_var*ri_x \n" ## Constrained to be >0
    )
    return(traitModel)
}


buildTraitY_l <- function(wavesList) {
    title <- "## Stable Trait/Random Intercepts for Y \n"
    ## First do indicator statements
    traitModel <- paste0(title, "ri_y =~ 1*ly", wavesList[1], " \n")
    for (w in wavesList[-1]) {
        traitModel <- paste0(
            traitModel,
            paste0("ri_y =~ 1*ly", w, " \n")
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
## Setup phantom variables for missing waves (if needed; returns comment if none)
################################################################################

buildPhantomX_l <- function(waves, xWaves) {
    if(length(xWaves) < length(c(1:waves))) {
        title <- "## Phantom X Variables\n"
        phantomX <- title
        phantomWaves <- c(1:waves)[-xWaves]
        for (w in phantomWaves) {
            phantomX <- paste0(phantomX, "lx", w, " =~ 0 \n")
        }
        return(phantomX)
    } else {
        return("## No Phantom X Variables;\n")
    }
}

buildPhantomY_l <- function(waves, yWaves) {
    if(length(yWaves) < length(c(1:waves))) {
        title <- "## Phantom Y Variables\n"
        phantomY <- title
        phantomWaves <- c(1:waves)[-yWaves]
        for (w in phantomWaves) {
            phantomY <- paste0(phantomY, "ly", w, " =~ 0\n")
        }
        return(phantomY)
    } else {
        return("## No Phantom Y Variables;\n")
    }
}



################################################################################
## Build latent occasion/observed part of model
################################################################################

buildObservedX_l <- function(xWaves, xIndicators) {
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

buildObservedY_l <- function(yWaves, yIndicators) {
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

buildStateX_l <- function(wavesList, stationarity = TRUE) {
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




buildStateY_l <- function(wavesList, stationarity = TRUE) {
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

buildLatentVarX_l <- function(waves) {
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

buildLatentVarY_l <- function(waves) {
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

buildARX_l <- function(waves, stationarity = TRUE) {
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
    if (stationarity == TRUE) {
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
                    " ~~ arv2x*arx",
                    w,
                    " \n"
                )
            )
        }
    } else {
        arXModel <- paste0(
            arXModel,
            subtitle2,
            "arx2 ~ arx1 \n"
        )
        for (w in wavesList[-c(1:2)]) {
            arXModel <- paste0(
                arXModel,
                paste0(
                    "arx",
                    w,
                    " ~ arx",
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
    }
    return(arXModel)
}

buildARY_l <- function(waves, stationarity = TRUE) {
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
    if (stationarity == TRUE) {
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
                    " ~~ arv2y*ary",
                    w,
                    " \n"
                )
            )
        }
    } else {
        arYModel <- paste0(
            arYModel,
            subtitle2,
            "ary2 ~ ary1 \n"
        )
        for (w in wavesList[-c(1:2)]) {
            arYModel <- paste0(
                arYModel,
                paste0(
                    "ary",
                    w,
                    " ~ ary",
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
    }
    return(arYModel)
}


################################################################################
## Build Correlations
################################################################################

buildTraitCors_l <- function() {
    traitCors <- "## Stable Trait Correlations \n ri_x ~~ cov_t*ri_y \n"
    return(traitCors)
}


buildArCors_l <- function(wavesList, stationarity = TRUE) {
    if (stationarity == TRUE) {
        arCors <- "## AR Correlations \narx1 ~~ cor_xy*ary1 \n"
        for (w in wavesList[-1]) {
            arCors <- paste0(
                arCors,
                paste0("arx", w, " ~~ cor_xyr*ary", w, " \n")
            )
        }
    } else {
        arCors <- "## AR Correlations \narx1 ~~ cor_xy*ary1 \n"
        for (w in wavesList[-1]) {
            arCors <- paste0(
                arCors,
                paste0("arx", w, " ~~ ary", w, " \n")
            )
        }
    }
    return(arCors)
}



buildStateCors_l <- function(wavesList, stationarity = TRUE) {
    if (stationarity == TRUE) {
        sCors <- "## State Correlations \n"
        for (w in wavesList) {
            sCors <- paste0(
                sCors,
                paste0("sx", w, " ~~ cor_s*sy", w, " \n")
            )
        }
    } else {
        sCors <- "## State Correlations \n"
        for (w in wavesList) {
            sCors <- paste0(
                sCors,
                paste0("sx", w, " ~~ sy", w, " \n")
            )
        }
    }
    return(sCors)
}

buildResidCorsX_l <- function(xWaves, xIndicators, constrainCors = TRUE) {
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


buildResidCorsY_l <- function(yWaves, yIndicators, constrainCors = TRUE) {
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

buildClYFromX_l <- function(wavesList, stationarity = TRUE) {
    title <- "### Cross-Lagged Path \n"
    subtitle <- "## Y predicted from X \n"
    clYonX <- paste0(title, subtitle)
    if (stationarity == TRUE) {
        for (w in wavesList[-1]) {
            clYonX <- paste0(
                clYonX, "ary", w, " ~ c*arx", (w - 1), " \n"
            )
        }
    } else {
        for (w in wavesList[-1]) {
            clYonX <- paste0(
                clYonX, "ary", w, " ~ arx", (w - 1), " \n"
            )
        }
    }
    return(clYonX)
}


buildClXFromY_l <- function(wavesList, stationarity = TRUE) {
    title <- "### Cross-Lagged Path \n"
    subtitle <- "## X predicted from Y \n"
    clXonY <- paste0(title, subtitle)
    if (stationarity == TRUE) {
        for (w in wavesList[-1]) {
            clXonY <- paste0(
                clXonY, "arx", w, " ~ d*ary", (w - 1), " \n"
            )
        }
    } else {
        for (w in wavesList[-1]) {
            clXonY <- paste0(
                clXonY, "arx", w, " ~ ary", (w - 1), " \n"
            )
        }
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

#' Build Lavaan Code for STARTS Model
#'
#' `buildLavaan()` produces lavaan code for variations of the general Stable Trait,
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
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
buildLavaan <- function(waves,
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
                        stationarity = TRUE,
                        constrainCors = TRUE,
                        limits = TRUE) {
    if (is.null(xWaves)) xWaves <- 1:waves
    if (is.null(yWaves)) yWaves <- 1:waves
    xWaves <- xWaves[xWaves <= waves]
    yWaves <- yWaves[yWaves <= waves]
    wavesList <- 1:waves
    startsModel <- "### Model Code Generated by buildModel.R"
    if (XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildPhantomX_l(waves, xWaves),
            buildObservedX_l(xWaves, xIndicators),
            buildLatentVarX_l(waves),
            sep=" \n")
    }
    if (YVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildPhantomY_l(waves, yWaves),
            buildObservedY_l(yWaves, yIndicators),
            buildLatentVarY_l(waves),
            sep = " \n"
        )
    }
    if (trait == TRUE & XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildTraitX_l(wavesList)
        )
    }
    if (trait == TRUE & YVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildTraitY_l(wavesList)
        )
    }
    if (state == TRUE & XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildStateX_l(wavesList, stationarity),
            sep = " \n"
        )
    }
    if (state == TRUE & YVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildStateY_l(wavesList, stationarity),
            sep = " \n"
        )
    }
    if (AR == TRUE & XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildARX_l(waves, stationarity),
            sep = " \n"
        )
    }
    if (AR == TRUE & YVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildARY_l(waves, stationarity),
            sep = " \n"
        )
    }
    if (crossLag == TRUE & YVar == TRUE  & XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildClYFromX_l(wavesList, stationarity),
            buildClXFromY_l(wavesList, stationarity),
            sep = " \n"
        )
    }
    if (trait == TRUE & YVar == TRUE & XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildTraitCors_l(),
            sep = " \n"
        )
    }
    if (AR == TRUE & YVar == TRUE & XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildArCors_l(wavesList, stationarity),
            sep = " \n"
        )
    }
    if (stateCor == TRUE &
        state == TRUE &
        YVar == TRUE &
        XVar == TRUE ) {
        startsModel <- paste(
            startsModel,
            buildStateCors_l(wavesList, stationarity),
            sep = " \n"
        )
    }

    if (XVar == TRUE &
        xIndicators > 1) {
        startsModel <- paste(
            startsModel,
            buildResidCorsX_l(xWaves, xIndicators, constrainCors),
            sep = " \n"
        )
    }

    if (YVar == TRUE &
        yIndicators > 1) {
        startsModel <- paste(
            startsModel,
            buildResidCorsY_l(yWaves, yIndicators, constrainCors),
            sep = " \n"
        )
    }

    
    ## Constraints
    cons_uni <- "arv2x == arvx - (a^2) * arvx \n"
    cons_y_uni <- "arv2y == arvy - (b^2) * arvy \n"
    cons_cl <- "arv2x == (1-a^2)*arvx - 2*a*d*cor_xy - d^2*arvy \n"
    cons_y_cl <- "arv2y == (1-b^2)*arvy - 2*b*c*cor_xy - c^2*arvx \n"
    cons_cor_noCl <- "cor_xyr == cor_xy - (a*cor_xy*b) \n"
    cons_cor_cl <- "cor_xyr == (1-a*b-c*d)*cor_xy - a*c*arvx - b*d*arvy \n"


    startsModel <- paste(
        startsModel,
        "## Constraints",
        sep = " \n"
    )

    ## Restrict variance to be non-zero
    if (limits == TRUE) {
        if (XVar == TRUE & trait == TRUE) {
            startsModel <- paste(
                startsModel,
                "x_t_var > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE & trait == TRUE) {
            startsModel <- paste(
                startsModel,
                "y_t_var > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE & AR == TRUE) {
            startsModel <- paste(
                startsModel,
                "arvx > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE &
            AR == TRUE &
            stationarity == TRUE) {
            startsModel <- paste(
                startsModel,
                "arv2x > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE & AR == TRUE) {
            startsModel <- paste(
                startsModel,
                "arvy > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE &
            AR == TRUE &
            stationarity == TRUE) {
            startsModel <- paste(
                startsModel,
                "arv2y > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE &
            state == TRUE &
            stationarity == TRUE) {
            startsModel <- paste(
                startsModel,
                "sx > 0",
                sep = " \n"
            )
        }
        if (YVar == TRUE &
            state == TRUE &
            stationarity == TRUE) {
            startsModel <- paste(
                startsModel,
                "sy > 0",
                sep = " \n"
            )
        }
        if (XVar == TRUE &
            YVar == TRUE) {
            if (trait == TRUE) {
                startsModel <- paste(
                    startsModel,
                    "cov_t < .99*(sqrt(x_t_var)*sqrt(y_t_var))",
                    "cov_t > -.99*(sqrt(x_t_var)*sqrt(y_t_var))",
                    sep = " \n"
                )
            }
            if (AR == TRUE) {
                startsModel <- paste(
                    startsModel,
                    "cor_xy < .99*(sqrt(arvx)*sqrt(arvy))",
                    "cor_xy > -.99*(sqrt(arvx)*sqrt(arvy))",
                    sep = " \n"
                )
            }
            if (state == TRUE & stateCor == TRUE) {
                startsModel <- paste(
                    startsModel,
                    "cor_s < .99*(sqrt(sx)*sqrt(sy))",
                    "cor_s > -.99*(sqrt(sx)*sqrt(sy))",
                    sep = " \n"
                )
            }
        }
    }
    if (stationarity == TRUE) {
        if (XVar == TRUE) {
            if (crossLag == FALSE | YVar == FALSE) {
                startsModel <- paste(
                    startsModel,
                    cons_uni,
                    sep = " \n"
                )
            } else if (crossLag == TRUE & YVar == TRUE) {
                startsModel <- paste(
                    startsModel,
                    cons_cl,
                    sep = " \n"
                )
            }
        }
        if (YVar == TRUE) {
            if (crossLag == FALSE | XVar == FALSE) {
                startsModel <- paste(
                    startsModel,
                    cons_y_uni,
                    sep = " \n"
                )
            } else if (crossLag == TRUE & XVar == TRUE) {
                startsModel <- paste(
                    startsModel,
                    cons_y_cl,
                    sep = " \n"
                )
            }
        }
        if (XVar == TRUE &
            YVar == TRUE &
            AR == TRUE ) {
            if (crossLag == TRUE) {
                startsModel <- paste(
                    startsModel,
                    cons_cor_cl,
                    sep = " \n"
                )
            } else {
                startsModel <- paste(
                    startsModel,
                    cons_cor_noCl,
                    sep = " \n"
                )
            }
        }
    }
    return(startsModel)
}


#' Builds Univariate Starts Model for X in Lavaan
#'
#' `lavaanStartsX()` Wrapper function for `run_starts_mplus()`. Produces lavaan
#' code for the univariate STARTS model (for variabes named X). 
#'
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @returns Returns character vector of code that can be run using lavaan.
#' @export
lavaanStartsX <- function(waves,
                          xWaves,
                          xIndicators = 1,
                          stationarity = TRUE,
                          constrainCors = TRUE,
                          limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        xIndicators = xIndicators,
        stationarity = stationarity,
        constrainCors = constrainCors,
        state = TRUE,
        YVar = FALSE,
        limits = limits
    )
}

#' Builds Univariate Starts Model for Y in Lavaan
#'
#' `lavaanStartsY()` Wrapper function for `run_starts_mplus()`. Produces lavaan
#' code for the univariate STARTS model (for variabes named Y). 
#'
#' @param waves Numeric value indicating total number of waves.
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @returns Returns character vector of code that can be run using lavaan.
#' @export
lavaanStartsY <- function(waves,
                          yWaves,
                          yIndicators = 1,
                          stationarity = TRUE,
                          constrainCors = TRUE,
                          limits = TRUE) {
    buildLavaan(
        waves = waves,
        yWaves = yWaves,
        yIndicators = yIndicators,
        stationarity = stationarity,
        constrainCors = constrainCors,
        state = TRUE,
        XVar = FALSE,
        limits = limits,
    )
}



#' Builds Bivariate Starts Model in Lavaan
#'
#' `lavaanStarts2()` Wrapper function for `run_starts_mplus()`. Produces lavaan
#' code for the standard bivariate STARTS model. 
#'
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param crossLag Logical value indicating whether to impose stationairty.
#'   Defaults to TRUE.
#' @param stateCor Logical value indicating whether to constrain correlations
#'   between wave-specific state variances to be equal across waves
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @returns Returns character vector of code that can be run using lavaan.
#' @export
lavaanStarts2 <- function(waves,
                          xWaves,
                          yWaves,
                          xIndicators = 1,
                          yIndicators = 1,
                          stationarity = TRUE,
                          constrainCors = TRUE,
                          crossLag = TRUE,
                          stateCor = TRUE,
                          limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        stationarity = stationarity,
        constrainCors = constrainCors,
        crossLag = crossLag,
        stateCor = stateCor,
        limits = limits
        )
}


#' Builds RI-CLPM Model in Lavaan
#'
#' `lavaanRiclpm()` Wrapper function for `run_starts_mplus()`. Produces lavaan
#' code for the standard Random-Intercept Cross-Lagged Panel Model. This is
#' equivalent to a STARTS model without the state variance.
#'
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @returns Returns character vector of code that can be run using lavaan.
#' @export
lavaanRiclpm <- function(waves,
                         xWaves = NULL,
                         yWaves = NULL,
                         xIndicators = xIndicators,
                         yIndicators = yIndicators,
                         stationarity = TRUE,
                         constrainCors = TRUE,
                         limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        stationarity = stationarity,
        constrainCors = constrainCors,
        state = FALSE,
        limits = limits
    )
}

#' Builds CLPM Model in Lavaan
#'
#' `lavaanjclpm()` Wrapper function for `run_starts_mplus()`. Produces lavaan
#' code for the standard cross-lagged panel model. This is equivalent to the
#' STARTS model with no stable trait or state variance.
#'
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @returns Returns character vector of code that can be run using lavaan.
#' @export
lavaanClpm <- function(waves,
                      xWaves = NULL,
                      yWaves = NULL,
                      xIndicators = 1,
                      yIndicators = 1,
                      stationarity = TRUE,
                      constrainCors = TRUE,
                      limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        trait = FALSE,
        state = FALSE,
        stateCor = FALSE,
        stationarity = stationarity,
        constrainCors = constrainCors,
        limits = limits
    )
}

#' Builds ARTS Model in Lavaan
#'
#' `lavaanArts()` Wrapper function for `run_starts_mplus()`. Produces lavaan
#' code for the autoregressive trait, state (ARTS) model. This is equivalent
#' to the STARTS model with no stable trait variance.
#'
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param stateCor Logical value indicating whether to include correlations
#'   between state components.
#' @param limits Logical value indicating whether to limit variances and
#'   correlations to valid values.
#' @returns Returns character vector of code that can be run using lavaan.
#' @export
lavaanArts <- function(waves,
                       xWaves = NULL,
                       yWaves = NULL,
                       xIndicators = 1,
                       yIndicators = 1,
                       stationarity = TRUE,
                       constrainCors = TRUE,
                       stateCor = FALSE,
                       limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        xIndicators = xIndicators,
        yIndicators = yIndicators,
        stationarity = stationarity,
        constrainCors = constrainCors,
        trait = FALSE,
        stateCor = stateCor,
        limits = limits
        )
}

summarizeLavaan <- function(modelStatement, fitObject) {
    ## Calculate values for summary
    parEst <- lavaan::parameterEstimates(fitObject)
    stdEst <- lavaan::standardizedSolution(fitObject)
    fit <- lavaan::fitMeasures(fitObject)
    ## Variance Decomp
    trait.x <- parEst[which(parEst$lhs=="ri_x" & parEst$rhs=="ri_x"), "est"]
    trait.y <- parEst[which(parEst$lhs=="ri_y" & parEst$rhs=="ri_y"), "est"]
    ar.x <- parEst[which(parEst$lhs=="arx1" & parEst$rhs=="arx1"), "est"]
    ar.y <- parEst[which(parEst$lhs=="ary1" & parEst$rhs=="ary1"), "est"]
    state.x <- parEst[which(parEst$lhs=="sx1" & parEst$rhs=="sx1"), "est"]
    state.y <- parEst[which(parEst$lhs=="sy1" & parEst$rhs=="sy1"), "est"]
    trait.x.p <- trait.x/(sum(trait.x, ar.x, state.x))
    trait.y.p <- trait.y/(sum(trait.y, ar.y, state.y))
    ar.x.p <- ar.x/(sum(trait.x, ar.x, state.x))
    ar.y.p <- ar.y/(sum(trait.y, ar.y, state.y))
    state.x.p <- state.x/(sum(trait.x, ar.x, state.x))
    state.y.p <- state.y/(sum(trait.y, ar.y, state.y))
    ## Correlations
    trait.cor <- stdEst[which(stdEst$lhs=="ri_x" & stdEst$rhs=="ri_y"), "est.std"]
    ar.cor <- stdEst[which(stdEst$lhs=="arx1" & stdEst$rhs=="ary1"), "est.std"]
    state.cor <- stdEst[which(stdEst$lhs=="sx1" & stdEst$rhs=="sy1"), "est.std"]
    ## Stability
    x.stab <- parEst[which(parEst$lhs=="arx2" & parEst$rhs=="arx1"), "est"]
    y.stab <- parEst[which(parEst$lhs=="ary2" & parEst$rhs=="ary1"), "est"]
    ## Cross-lags
    yOnX <- parEst[which(parEst$lhs=="ary2" & parEst$rhs=="arx1"), "est"]
    xOnY <- parEst[which(parEst$lhs=="arx2" & parEst$rhs=="ary1"), "est"]
    ## Fit
    chi2 <- fit[["chisq"]]
    chi2df <- fit[["df"]]
    chi2p <- fit[["pvalue"]]
    cfi <- fit[["cfi"]]
    tli <- fit[["tli"]]
    srmr <- fit[["srmr"]]
    rmsea <- fit[["rmsea"]]
    outputList <- list(
        modelStatement = modelStatement,
        lavaanFit = fitObject,
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
    class(outputList) <- "pclObject"
    return(outputList)
}




#' Build and Run Lavaan Code for STARTS Model
#'
#' `run_starts_lavaan()` produces and runs lavaan code for variations of the
#' general Stable Trait, Autoregressive Trait, State Model. Univariate and
#' bivariate versions can be specified. Reduced versions of the model can be
#' specified by setting variance components for trait, AR, or state to be
#' `FALSE`.
#'
#' @param data Dataframe with all variables.
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
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @param ... Additional arguments to the `lavaan()` command.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
run_starts_lavaan <- function(data,
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
                              stateCor = FALSE,
                              stationarity = TRUE,
                              constrainCors = TRUE,
                              limits = TRUE,
                              ...) {
    startsModel <- buildLavaan(waves = waves,
                               XVar = XVar,
                               YVar = YVar,
                               xWaves = xWaves,
                               yWaves = yWaves,
                               xIndicators = xIndicators,
                               yIndicators = yIndicators,
                               trait = trait,
                               AR = AR,
                               state = state,
                               crossLag = crossLag,
                               stateCor = stateCor,
                               stationarity = stationarity,
                               constrainCors = constrainCors,
                               limits = limits)
    startsFit <- lavaan::lavaan(startsModel, data, ...)
    outputList <- summarizeLavaan(startsModel, startsFit)
    summary(outputList)
    return(outputList)
}


#' Build and Run Lavaan Code for RI-CLPM Model
#'
#' `run_riclpm_lavaan()` produces and runs lavaan code for the random-intercept
#' cross-lagged panel model. 
#'
#' @param data Dataframe with all variables.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @param ... Additional arguments to the `lavaan()` command.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
run_riclpm_lavaan <- function(data,
                              waves,
                              xWaves = NULL,
                              yWaves = NULL,
                              xIndicators = 1,
                              yIndicators = 1,
                              stationarity = TRUE,
                              constrainCors = TRUE,
                              limits = TRUE,
                              ...) {
    startsModel <- lavaanRiclpm(waves = waves,
                                xWaves = xWaves,
                                yWaves = yWaves,
                                xIndicators = xIndicators,
                                yIndicators = yIndicators,
                                stationarity = stationarity,
                                constrainCors = constrainCors,
                                limits = limits)
    startsFit <- lavaan::lavaan(startsModel, data, ...)
    outputList <- summarizeLavaan(startsModel, startsFit)
    summary(outputList)
    return(outputList)
}

#' Build and Run Lavaan Code for CLPM Model
#'
#' `run_clpm_lavaan()` produces and runs lavaan code for the standard lag-1
#' cross-lagged panel model.
#'
#' @param data Dataframe with all variables.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @param ... Additional arguments to the `lavaan()` command.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
run_clpm_lavaan <- function(data,
                              waves,
                              xWaves = NULL,
                              yWaves = NULL,
                              xIndicators = 1,
                              yIndicators = 1,
                              stationarity = TRUE,
                              constrainCors = TRUE,
                              limits = TRUE,
                              ...) {
    startsModel <- lavaanClpm(waves = waves,
                              xWaves = xWaves,
                              yWaves = yWaves,
                              xIndicators = xIndicators,
                              yIndicators = yIndicators,
                              stationarity = stationarity,
                              constrainCors = constrainCors,
                              limits = limits)
    startsFit <- lavaan::lavaan(startsModel, data, ...)
    outputList <- summarizeLavaan(startsModel, startsFit)
    summary(outputList)
    return(outputList)
}

#' Build and Run Lavaan Code for ARTS Model
#'
#' `run_arts_lavaan()` produces and runs lavaan code for the ARTS model, which
#' includes autoregressive and state variance, but no trait variance. This is
#' equivalent to a 'factor CLPM'.
#'
#' @param data Dataframe with all variables.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1. 
#' @param stateCor logical value indicating whether to include state
#'   correlations. Defaults to TRUE.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @param ... Additional arguments to the `lavaan()` command.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
run_arts_lavaan <- function(data,
                            waves,
                            xWaves = NULL,
                            yWaves = NULL,
                            xIndicators = 1,
                            yIndicators = 1,
                            stateCor = FALSE,
                            stationarity = TRUE,
                            constrainCors = TRUE,
                            limits = TRUE,
                            ...) {
    startsModel <- lavaanArts(waves = waves,
                              xWaves = xWaves,
                              yWaves = yWaves,
                              xIndicators = xIndicators,
                              yIndicators = yIndicators,
                              stateCor = stateCor,
                              stationarity = stationarity,
                              constrainCors = constrainCors,
                              limits = limits)
    startsFit <- lavaan::lavaan(startsModel, data, ...)
    outputList <- summarizeLavaan(startsModel, startsFit)
    summary(outputList)
    return(outputList)
}


#' Build and Run Lavaan Code for Univariate STARTS Model
#'
#' `run_startsX_lavaan()` produces and runs lavaan code for the univariate
#' STARTS models (for X variables). 
#'
#' @param data Dataframe with all variables.
#' @param waves Numeric value indicating total number of waves.
#' @param xWaves Vector of actual waves for X (omit if same as waves).
#' @param xIndicators Numeric value indicatoring number of indicators for X.
#'   Defaults to 1.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @param ... Additional arguments to the `lavaan()` command.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
run_startsX_lavaan <- function(data,
                              waves,
                              xWaves = NULL,
                              xIndicators = 1,
                              stationarity = TRUE,
                              constrainCors = TRUE,
                              limits = TRUE,
                              ...) {
    startsModel <- lavaanStartsX(waves = waves,
                                 xWaves = xWaves,
                                 xIndicators = xIndicators,
                                 stationarity = stationarity,
                                 constrainCors = constrainCors,
                                 limits = limits)
    startsFit <- lavaan::lavaan(startsModel, data, ...)
    outputList <- summarizeLavaan(startsModel, startsFit)
    summary(outputList)
    return(outputList)
}


#' Build and Run Lavaan Code for Univariate STARTS Model
#'
#' `run_startsY_lavaan()` produces and runs lavaan code for the univariate
#' STARTS models (for Y variables). 
#'
#' @param data Dataframe with all variables.
#' @param waves Numeric value indicating total number of waves.
#' @param yWaves Vector of actual waves for Y (omit if same as waves).
#' @param yIndicators Numeric value indicatoring number of indicators for Y.
#'   Defaults to 1.
#' @param stationarity logical value indicating whether to impose stationarity.
#'   Defaults to TRUE.
#' @param constrainCors logical value indicating whether to constrain
#'   correlations between same indicator at different waves (when there are more
#'   than one indicator.
#' @param limits Logical value indicating whether to constrain variances and
#'   correlations to possible values.
#' @param ... Additional arguments to the `lavaan()` command.
#' @returns Character vector representing the Mplus code for the model statement.
#' @export
run_startsY_lavaan <- function(data,
                              waves,
                              yWaves = NULL,
                              yIndicators = 1,
                              stationarity = TRUE,
                              constrainCors = TRUE,
                              limits = TRUE,
                              ...) {
    startsModel <- lavaanStartsY(waves = waves,
                                 yWaves = yWaves,
                                 yIndicators = yIndicators,
                                 stationarity = stationarity,
                                 constrainCors = constrainCors,
                                 limits = limits)
    startsFit <- lavaan::lavaan(startsModel, data, ...)
    outputList <- summarizeLavaan(startsModel, startsFit)
    summary(outputList)
    return(outputList)
}



#' Summarizes results from `run_starts_lavaan()`
#'
#' @param object Results from `run_starts_lavaan()`.
#' @param ... Additional arguments to `summary()`.
#'
#' @export
summary.pclObject <- function(object, ...) {
    cat("Model Summary: \n")
    cat("Model: Chi2 (df = ",
        sprintf("%i", object$chi2df), ") = ",
        sprintf("%.3f", object$chi2), ", p = ",
        sprintf("%.3f", object$chi2p), "\n")
    cat("\n")
    cat("Fit Indices:")
    cat("\n")
    cat("CFI = ",
        sprintf("%.3f", object$cfi),
        ", TLI = ",
        sprintf("%.3f", object$tli),
        ", SRMR = ",
        sprintf("%.3f", object$srmr), "\n")
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
        "\n")
    cat("Y Variable: \n")
    cat("Trait: ",
        sprintf("%.3f", object$trait.y),
        ", Autoregressive: ",
        sprintf("%.3f", object$ar.y),
        ", State: ",
        sprintf("%.3f", object$state.y), "\n")
    cat("\n")
    cat("Stability: \n")
    cat("X: ", sprintf("%.3f", object$x.stab), "\n")
    cat("Y: ", sprintf("%.3f", object$y.stab), "\n")
    cat("\n")
    cat("Cross-Lag Paths: \n")
    cat("Y predicted from X: ", sprintf("%.3f", object$yOnX) ,"\n")
    cat("X predicted from Y: ", sprintf("%.3f", object$xOnY),"\n")
    cat("Correlations: \n")
    cat("\n")
    cat("Stable Trait: ", sprintf("%.3f", object$trait.cor), "\n")
    cat("Autoregressive Trait: ", sprintf("%.3f", object$ar.cor), "\n")
    cat("State: ", sprintf("%.3f", object$state.cor), "\n")
    }
