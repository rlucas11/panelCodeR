compareUnivariate <- function(data, waves, xWaves=NULL) {
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
        if (!is.null(fit.est[[which(fit.est$param == "TRAIT.P."), "est"]])) {
            singleModelResults[1, 1] <- fit.est[[which(fit.est$param == "TRAIT.P."), "est"]]
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
        if (!is.null(fit$results$summaries$RMSEA_Estimate)) {
            singleModelResults[1, 9] <- fit$results$summaries$RMSEA_Estimate
        } else {
            singleModelResults[1, 9] <- NA
        }
        return(singleModelResults)
    }
    if (is.null(xWaves)) xWaves <- 1:waves
    starts <- run_startsx_mplus(
        data,
        waves,
        xWaves,
        title="starts"
    )
    results[1, 1:9] <- unlist(collectResults(starts)[1, 1:9])
    st <- run_starts_mplus(
        data,
        waves,
        xWaves,
        AR = FALSE,
        YVar = FALSE,
        title = "st"
    )
    results[1, 10:18] <- unlist(collectResults(st)[1, 1:9])
    arts <- run_starts_mplus(
        data,
        waves,
        xWaves,
        trait = FALSE,
        YVar = FALSE,
        title = "arts"
    )
    results[1, 19:27] <- unlist(collectResults(arts)[1, 1:9])
    start <- run_starts_mplus(
        data,
        waves,
        xWaves,
        YVar = FALSE,
        state = FALSE,
        title = "start"
    )
    results[1, 28:36] <- unlist(collectResults(start)[1, 1:9])
    art <- run_starts_mplus(
        data,
        waves,
        xWaves,
        trait = FALSE,
        YVar = FALSE,
        state = FALSE,
        title = "art"
    )
    results[1, 37:45] <- unlist(collectResults(art)[1, 1:9])
    aicCols <- paste("aic", c("starts", "st", "arts", "start", "art"), sep = ".")
    bicCols <- paste("bic", c("starts", "st", "arts", "start", "art"), sep = ".")
    results$best.aic <- colnames(results[, aicCols])[apply(results[, aicCols], 1, which.min)]
    results$best.bic <- colnames(results[, bicCols])[apply(results[, bicCols], 1, which.min)]
    return(results)
}
