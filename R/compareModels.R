#' Compares Univariate Variations of STARTS Model (Old Version)
#'
#' `compareUnivariate()` produces code, runs, and compares variations of the
#' Stable Trait, Autoregressive Trait, State Model. Specifically, it
#' sequentially drops each component (stable trait, autoregressive trait, and
#' state) and then uses BIC to identify the best model. 
#'
#' @param data Dataframe with multiwave data. Variables should be names 'x1' to
#'   'xw', where 'w' is the number of waves.
#' @param waves Numeric value indicating the number of waves.
#' @param xWaves Numeric vector indicating which waves actually exist (e.g.,
#'   'c(1:3, 5:7)' if there are 7 waves but Wave 6 is missing.
#' @returns Matrix of results for each model. This includes basic fit indices
#'   as well as variance decomposition and stability estimates. 
#' @export
compareUnivariate_old <- function(data, waves, xWaves=NULL) {
    results <- data.frame(
        trait.p.starts = numeric(),
        ar.p.starts = numeric(),
        st.p.starts = numeric(),
        stability.starts = numeric(),
        aic.starts = numeric(),
        bic.starts = numeric(),
        chisq.starts = numeric(),
        df.starts = numeric(),
        rmsea.starts = numeric(),
        trait.p.st = numeric(),
        ar.p.st = numeric(),
        st.p.st = numeric(),
        stability.st = numeric(),
        aic.st = numeric(),
        bic.st = numeric(),
        chisq.st = numeric(),
        df.st = numeric(),
        rmsea.st = numeric(),
        trait.p.arts = numeric(),
        ar.p.arts = numeric(),
        st.p.arts = numeric(),
        stability.arts = numeric(),
        aic.arts = numeric(),
        bic.arts = numeric(),
        chisq.arts = numeric(),
        df.arts = numeric(),
        rmsea.arts = numeric(),
        trait.p.start = numeric(),
        ar.p.start = numeric(),
        st.p.start = numeric(),
        stability.start = numeric(),
        aic.start = numeric(),
        bic.start = numeric(),
        chisq.start = numeric(),
        df.start = numeric(),
        rmsea.start = numeric(),
        trait.p.art = numeric(),
        ar.p.art = numeric(),
        st.p.art = numeric(),
        stability.art = numeric(),
        aic.art = numeric(),
        bic.art = numeric(),
        chisq.art = numeric(),
        df.art = numeric(),
        rmsea.art = numeric(),
        best.aic = character(),
        best.bic = character()
    )
    collectResults <- function(fit) {
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
            rmsea = numeric()
        )
        fit.est <- fit$results$parameters$unstandardized
        if (!is.null(fit.est[[which(fit.est$param == "TRAIT.X"), "est"]])) {
            singleModelResults[1, 1] <- fit.est[[which(fit.est$param == "TRAIT.X"), "est"]]
        } else {
            singleModelResults[1, 1] <- NA
        }
        if (!is.null(fit.est[[which(fit.est$param == "AR.P.X"), "est"]])) {
            singleModelResults[1, 2] <- fit.est[[which(fit.est$param == "AR.P.X"), "est"]]
        } else {
            singleModelResults[1, 2] <- NA
        }
        if (!is.null(fit.est[[which(fit.est$param == "ST.P.X"), "est"]])) {
            singleModelResults[1, 3] <- fit.est[[which(fit.est$param == "ST.P.X"), "est"]]
        } else {
            singleModelResults[1, 3] <- NA
        }
        if (!is.null(fit.est[[which(fit.est$param == "STAB.X"), "est"]])) {
            singleModelResults[1, 4] <- fit.est[[which(fit.est$param == "STAB.X"), "est"]]
        } else {
            singleModelResults[1, 4] <- NA
        }
        if (!is.null(fit$results$summaries$AIC)) {
            singleModelResults[1, 5] <- fit$results$summaries$AIC
        } else {
            singleModelResults[1, 5] <- NA
        }
        if (!is.null(fit$results$summaries$BIC)) {
            singleModelResults[1, 6] <- fit$results$summaries$BIC
        } else {
            singleModelResults[1, 6] <- NA
        }
        if (!is.null(fit$results$summaries$ChiSqM_Value)) {
            singleModelResults[1, 7] <- fit$results$summaries$ChiSqM_Value
        } else {
            singleModelResults[1, 7] <- NA
        }
        if (!is.null(fit$results$summaries$ChiSqM_DF)) {
            singleModelResults[1, 8] <- fit$results$summaries$ChiSqM_DF
        } else {
            singleModelResults[1, 8] <- NA
        }
        if (!is.null(fit$results$summaries$ChiSqM_Value)) {
            singleModelResults[1, 7] <- fit$results$summaries$ChiSqM_Value
        } else {
            singleModelResults[1, 7] <- NA
        }
        if (!is.null(fit$results$summaries$CFI)) {
            singleModelResults[1, 9] <- fit$results$summaries$CFI
        } else {
            singleModelResults[1, 9] <- NA
        }
        if (!is.null(fit$results$summaries$RMSEA_Estimate)) {
            singleModelResults[1, 10] <- fit$results$summaries$RMSEA_Estimate
        } else {
            singleModelResults[1, 10] <- NA
        }
        return(singleModelResults)
    }
    if (is.null(xWaves)) xWaves <- 1:waves
    starts <- run_startsx_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        title="starts"
    )
    results[1, 1:9] <- unlist(collectResults(starts$mplusOutput)[1, 1:9])
    st <- run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        AR = FALSE,
        YVar = FALSE,
        title = "st"
    )
    results[1, 10:18] <- unlist(collectResults(st$mplusOutput)[1, 1:9])
    arts <- run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        trait = FALSE,
        YVar = FALSE,
        title = "arts"
    )
    results[1, 19:27] <- unlist(collectResults(arts$mplusOutput)[1, 1:9])
    start <- run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        YVar = FALSE,
        state = FALSE,
        title = "start"
    )
    results[1, 28:36] <- unlist(collectResults(start$mplusOutput)[1, 1:9])
    art <- run_starts_mplus(
        data = data,
        waves = waves,
        xWaves = xWaves,
        trait = FALSE,
        YVar = FALSE,
        state = FALSE,
        title = "art"
    )
    results[1, 37:45] <- unlist(collectResults(art$mplusOutput)[1, 1:9])
    aicCols <- paste("aic", c("starts", "st", "arts", "start", "art"), sep = ".")
    bicCols <- paste("bic", c("starts", "st", "arts", "start", "art"), sep = ".")
    results$best.aic <- colnames(results[, aicCols])[apply(results[, aicCols], 1, which.min)]
    results$best.bic <- colnames(results[, bicCols])[apply(results[, bicCols], 1, which.min)]
    return(results)
}

#' Compares Univariate Variations of STARTS Model
#'
#' `compareUnivariate()` produces code, runs, and compares variations of the
#' Stable Trait, Autoregressive Trait, State Model. Specifically, it
#' sequentially drops each component (stable trait, autoregressive trait, and
#' state) and then uses BIC to identify the best model. 
#'
#' @param data Dataframe with multiwave data. Variables should be names 'x1' to
#'   'xw', where 'w' is the number of waves.
#' @param program Name of the program. Defaults to mplus as it is much quicker.
#' @param title Optional character string to be attended to title of model.
#' @export
compareUnivariate <- function(data, program = "mplus", title=NULL) {
    results <- data.frame(
        trait.p.starts = numeric(),
        ar.p.starts = numeric(),
        st.p.starts = numeric(),
        stability.starts = numeric(),
        aic.starts = numeric(),
        bic.starts = numeric(),
        chisq.starts = numeric(),
        df.starts = numeric(),
        cfi.starts = numeric(),
        tli.starts = numeric(),
        srmr.starts = numeric(),
        rmsea.starts = numeric(),
        trait.p.st = numeric(),
        ar.p.st = numeric(),
        st.p.st = numeric(),
        stability.st = numeric(),
        aic.st = numeric(),
        bic.st = numeric(),
        chisq.st = numeric(),
        df.st = numeric(),
        cfi.st = numeric(),
        tli.st = numeric(),
        srmr.st = numeric(),
        rmsea.st = numeric(),
        trait.p.arts = numeric(),
        ar.p.arts = numeric(),
        st.p.arts = numeric(),
        stability.arts = numeric(),
        aic.arts = numeric(),
        bic.arts = numeric(),
        chisq.arts = numeric(),
        df.arts = numeric(),
        cfi.arts = numeric(),
        tli.arts = numeric(),
        srmr.arts = numeric(),
        rmsea.arts = numeric(),
        trait.p.start = numeric(),
        ar.p.start = numeric(),
        st.p.start = numeric(),
        stability.start = numeric(),
        aic.start = numeric(),
        bic.start = numeric(),
        chisq.start = numeric(),
        df.start = numeric(),
        cfi.start = numeric(),
        tli.start = numeric(),
        srmr.start = numeric(),
        rmsea.start = numeric(),
        trait.p.art = numeric(),
        ar.p.art = numeric(),
        st.p.art = numeric(),
        stability.art = numeric(),
        aic.art = numeric(),
        bic.art = numeric(),
        chisq.art = numeric(),
        df.art = numeric(),
        cfi.art = numeric(),
        tli.art = numeric(),
        srmr.art = numeric(),
        rmsea.art = numeric(),
        best.aic = character(),
        best.bic = character()
    )
    
    starts <- tryCatch(
        panelcoder(
            data = data,
            title = paste(title, "starts", sep = "_"),
            program = program
        ),
        error = function(e) {
            return(NA)
        },
        warning = function(w) {
            return(NA)
        }
    )
    if (!is.na(starts[1])) {
        results[1, 1:12] <- unlist(.collectResults(starts)[1, 1:12])
    } else {
        results[1, 1:12] <- NA
    }
    
    st <- tryCatch(
        panelcoder(
            data = data,
            panelModel = "sts",
            program = program,
            title = paste(title, "sts", sep = "_")
        ),
        error = function(e) {
            return(NA)
        },
        warning = function(w) {
            return(NA)
        }
    )
    if (!is.na(st[1])) {
        results[1, 13:24] <- unlist(.collectResults(st)[1, 1:12])
    } else {
        results[1, 13:24] <- NA
    }
    arts <- tryCatch(
        panelcoder(
            data = data,
            panelModel = "arts",
            program = program,
            title = paste(title, "arts", sep = "_")
        ),
        error = function(e) {
            return(NA)
        },
        warning = function(w) {
            return(NA)
        }
    )
    if (!is.na(arts[1])) {
        results[1, 25:36] <- unlist(.collectResults(arts)[1, 1:12])
    } else {
        results[1, 25:36] <- NA
    }
    start <- tryCatch(panelcoder(
            data = data,
            panelModel = "riclpm",
            program = program,
            title = paste(title, "start", sep = "_")
        ),
        error = function(e) {
            return(NA)
        },
        warning = function(w) {
            return(NA)
        }
    )
    if (!is.na(start[1])) {
        results[1, 37:48] <- unlist(.collectResults(start)[1, 1:12])
    } else {
        results[1, 37:48] <- NA
    }
        art <- tryCatch(
            panelcoder(
                data = data,
                panelModel = "clpm",
                program = program,
                title = paste(title, "art", sep = "_")
            ),
            error = function(e) {
                return(NA)
            },
            warning = function(w) {
                return(NA)
            }
        )
    if (!is.na(art[1])) {
        results[1, 49:60] <- unlist(.collectResults(art)[1, 1:12])
    } else {
        results[1, 49:60] <- NA
    }
    
    aicCols <- paste("aic", c("starts", "st", "arts", "start", "art"), sep = ".")
    bicCols <- paste("bic", c("starts", "st", "arts", "start", "art"), sep = ".")
    results$best.aic <- colnames(results[, aicCols])[apply(results[, aicCols], 1, which.min)]
    results$best.bic <- colnames(results[, bicCols])[apply(results[, bicCols], 1, which.min)]
    return(results)
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
    
