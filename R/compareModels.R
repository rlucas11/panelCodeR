#' Compares Univariate Variations of Models Estimated Using panelCodeR
#'
#' `compareUnivariate()` produces code for, runs, and compares various models for
#' analyzing panel data. It then uses BIC to identify the best model. 
#'
#' @param data Dataframe with multiwave data. 
#' @param program Name of the program. Defaults to mplus as it is much quicker.
#' @param models List of models to compare. Defaults to list("starts", "riclpm",
#'   "arts", "sts", "clpm"). Can also include "dpm_c", "dpm_p", "gclm", "alt",
#'   "lgm_sm", or "lgcm" (see panelCodeR documentation). 
#' @param title Optional character string to be attended to title of model.
#' @export
compareUnivariate <- function(data,
                              program = "mplus",
                              models = list("starts",
                                            "riclpm",
                                            "arts",
                                            "sts",
                                            "clpm"),
                              title=NULL) {

    compareOut <- list()
    for (i in 1:length(models)) {
        pcOut <- tryCatch(
            panelcoder(
                data = data,
                title = paste(title, models[[i]], sep = "_"),
                panelModel = models[[i]],
                program = program
            ),
            error = function(e) {
                return(NA)
            },
            warning = function(w) {
                return(NA)
            }
        )
        if (!is.na(pcOut[1])) {
            compareOut[i] <- list(.collectResults(pcOut)[1, 1:12])
        } else {
            compareOut[i] <- NA
        }
    }
    return(compareOut)
}


    
    


.collectResults <- function(fit) {
    singleModelResults <- data.frame(
        trait = numeric(),
        ar = numeric(),
        st = numeric(),
        stab = numeric(),
        aic = numeric(),
        bic = numeric(),
        chisq = numeric(),
        df = numeric(),
        cfi = numeric(),
        tli = numeric(),
        srmr = numeric(),
        rmsea = numeric()
    )
    pcResults <- fit[[1]]
    if (!is.null(pcResults$trait.x)) {
        singleModelResults[1, 1] <- pcResults$trait.x
    } else {
        singleModelResults[1, 1] <- NA
    }
    if (!is.null(pcResults$ar.x)) {
        singleModelResults[1, 2] <- pcResults$ar.x
    } else {
        singleModelResults[1, 2] <- NA
    }
    if (!is.null(pcResults$state.x)) {
        singleModelResults[1, 3] <- pcResults$state.x
    } else {
        singleModelResults[1, 3] <- NA
    }
    if (!is.null(pcResults$x.stab)) {
        singleModelResults[1, 4] <- pcResults$x.stab
    } else {
        singleModelResults[1, 4] <- NA
    }
    if (!is.null(pcResults$aic)) {
        singleModelResults[1, 5] <- pcResults$aic
    } else {
        singleModelResults[1, 5] <- NA
    }
    if (!is.null(pcResults$bic)) {
        singleModelResults[1, 6] <- pcResults$bic
    } else {
        singleModelResults[1, 6] <- NA
    }
    if (!is.null(pcResults$chi2)) {
        singleModelResults[1, 7] <- pcResults$chi2
    } else {
        singleModelResults[1, 7] <- NA
    }
    if (!is.null(pcResults$chi2df)) {
        singleModelResults[1, 8] <- pcResults$chi2df
    } else {
        singleModelResults[1, 8] <- NA
    }
    if (!is.null(pcResults$chi2p)) {
        singleModelResults[1, 7] <- pcResults$chi2p
    } else {
        singleModelResults[1, 7] <- NA
    }
    if (!is.null(pcResults$cfi)) {
        singleModelResults[1, 9] <- pcResults$cfi
    } else {
        singleModelResults[1, 9] <- NA
    }
    if (!is.null(pcResults$tli)) {
        singleModelResults[1, 10] <- pcResults$tli
    } else {
        singleModelResults[1, 10] <- NA
    }
    if (!is.null(pcResults$srmr)) {
        singleModelResults[1, 11] <- pcResults$srmr
    } else {
        singleModelResults[1, 11] <- NA
    }
    if (!is.null(pcResults$rmsea)) {
        singleModelResults[1, 12] <- pcResults$rmsea
    } else {
        singleModelResults[1, 12] <- NA
    }
    return(singleModelResults)
}
    
