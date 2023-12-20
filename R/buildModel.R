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

buildObservedX_l <- function(waves, xWaves) {
    title <- "### Latent Occasion Variables for X \n"
    ## First do indicator statements
    observedModel <- title
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
        "## Latent Occasion Variance \n"
    )
    for (w in 1:waves) {
        observedModel <- paste0(
            observedModel,
            "lx",
            w,
            " ~~ ",
            0,
            "*lx",
            w,
            " \n"
        )
    }
    return(observedModel)
}


buildObservedY_l <- function(waves, yWaves) {
    title <- "### Latent Occasions Variables for Y \n"
    observedModel <- title
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
        "## Latent Occasion Variance \n"
    )
    for (w in 1:waves) {
        observedModel <- paste0(
            observedModel,
            "ly",
            w,
            " ~~ ",
            0,
            "*ly",
            w,
            " \n"
        )
    }
    return(observedModel)
}

################################################################################
## Build State Variance
################################################################################

buildStateX_l <- function(wavesList, stationarity = TRUE) {
    title <- "\n### X States;\n"
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
        "\n## State Variance;\n"
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
    title <- "\n### Y States;\n"
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
        "\n## State Variance;\n"
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
## Build autoregressive part of model
################################################################################

buildARX_l <- function(waves, stationarity = TRUE) {
    wavesList <- 1:waves
    title <- "### Autoregressive Part for X;\n"
    subtitle1 <- "## Indicator Statements;\n"
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
    title <- "### Autoregressive Part for Y;\n"
    subtitle1 <- "## Indicator Statements;\n"
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

buildLavaan <- function(waves,
                       xWaves = NULL,
                       yWaves = NULL,
                       stationarity = TRUE,
                       trait = TRUE,
                       XVar = TRUE,
                       YVar = TRUE,
                       crossLag = TRUE,
                       state = TRUE,
                       AR = TRUE,
                       stateCor = FALSE,
                       limits = TRUE) {
    if (is.null(xWaves)) xWaves <- 1:waves
    if (is.null(yWaves)) yWaves <- 1:waves
    wavesList <- 1:waves
    startsModel <- "### Model Code Generated by buildModel.R"
    if (XVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildPhantomX_l(waves, xWaves),
            buildObservedX_l(waves, xWaves),
            sep=" \n")
    }
    if (YVar == TRUE) {
        startsModel <- paste(
            startsModel,
            buildPhantomY_l(waves, yWaves),
            buildObservedY_l(waves, yWaves),
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


lavaanStartsX <- function(waves,
                         xWaves,
                         stationarity = TRUE,
                         limits = TRUE) {
    buildLavaan(
        waves,
        xWaves,
        stationarity = stationarity,
        state = TRUE,
        YVar = FALSE,
        limits = limits
    )
}

lavaanStartsY <- function(waves,
                         yWaves,
                         stationarity = TRUE,
                         limits = TRUE) {
    buildLavaan(
        waves,
        yWaves,
        stationarity = stationarity,
        state = TRUE,
        XVar = FALSE,
        limits = limits,
    )
}



lavaanStarts2 <- function(waves,
                         xWaves,
                         yWaves,
                         stationarity = TRUE,
                         crossLag = TRUE,
                         stateCor = TRUE,
                         limits = TRUE) {
    buildLavaan(
        waves,
        xWaves,
        yWaves,
        stationarity = stationarity,
        state = TRUE,
        crossLag = crossLag,
        stateCor = stateCor,
        limits = limits
        )
}


lavaanRiclpm <- function(waves,
                        xWaves = NULL,
                        yWaves = NULL,
                        limits = TRUE,
                        stationarity = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        stationarity = stationarity,
        state = FALSE,
        limits = limits
    )
}

lavaanClpm <- function(waves,
                      xWaves = NULL,
                      yWaves = NULL,
                      stationarity = TRUE,
                      limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        trait = FALSE,
        state = FALSE,
        stationarity = stationarity,
        limits = limits
    )
}

lavaanArts2 <- function(waves,
                       xWaves = NULL,
                       yWaves = NULL,
                       stationarity = TRUE,
                       crossLag = TRUE,
                       stateCor = TRUE,
                       limits = TRUE) {
    buildLavaan(
        waves = waves,
        xWaves = xWaves,
        yWaves = yWaves,
        stationarity = stationarity,
        trait = FALSE,
        crossLag = crossLag,
        stateCor = stateCor,
        limits = limits
        )
}
