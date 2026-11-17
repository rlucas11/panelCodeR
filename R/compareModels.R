#' Run and Compare Models Estimated Using panelCodeR
#'
#' `compareModels()` produces code for, runs, and compares various models for
#' analyzing panel data. It produces a table of basic results. 
#'
#' @param data Dataframe with multiwave data. 
#' @param program Name of the program. Defaults to mplus as it is much quicker.
#' @param models List of models to compare. Defaults to list("starts", "riclpm",
#'   "arts", "sts", "clpm"). Can also include "dpm_c", "dpm_p", "gclm", "alt",
#'   "lcmsr", or "lgcm" (see panelCodeR documentation). 
#' @param title Optional character string to be attended to title of model.
#' @param ... Options to pass to panelcoder function.
#' @export
compareModels <- function(data,
                          program = "mplus",
                          models = list("starts",
                                        "riclpm",
                                        "arts",
                                        "sts",
                                        "clpm"),
                          title=NULL,
                          ...) {

    compareOut <- list()
    for (i in 1:length(models)) {
        pcOut <- tryCatch(
            panelcoder(
                data = data,
                title = paste(title, models[[i]], sep = "_"),
                panelModel = models[[i]],
                program = program,
                ...
            ),
            error = function(e) {
                return(NA)
            },
            warning = function(w) {
                return(NA)
            }
        )
        if (!is.na(pcOut[1])) {
            compareOut[i] <- list(.collectResults(pcOut))
        } else {
            compareOut[i] <- list(.collectResults(list(NULL)))
        }
        compareOut[[i]]$model <- models[[i]]
    }
    compareOut <- do.call(rbind, compareOut)
    print(compareOut)
    return(compareOut)
}


    
    


.collectResults <- function(fit) {
    fields <- c(trait = "trait.x", ar = "ar.x", st = "state.x",
                stab = "x.stab", aic = "aic", bic = "bic",
                chisq = "chi2", df = "chi2df", chisqp = "chi2p",
                cfi = "cfi", tli = "tli", srmr = "srmr", rmsea = "rmsea")
    pcResults <- fit[[1]]
    values <- lapply(fields, function(f) {
        v <- if (!is.null(pcResults)) pcResults[[f]]
        if (is.null(v)) NA_real_ else v
    })
    out <- as.data.frame(values)
    out$warning <- if (length(fit) >= 7 && !is.null(fit[[7]])) fit[[7]] else NA_character_
    out$error <- if (length(fit) >= 8 && !is.null(fit[[8]])) fit[[8]] else NA_character_

    ## Diagnostics (NA when the result has none, e.g., a failed run)
    status <- if (is.list(fit)) fit[["status"]] else NULL
    diag <- if (is.list(fit)) fit[["diagnostics"]] else NULL
    flag <- function(name, na) {
        v <- if (!is.null(status)) status[[name]]
        if (is.null(v) || length(v) == 0) na else v[1]
    }
    out$converged <- flag("converged", NA)
    out$seTrusted <- flag("seTrusted", NA)
    out$admissible <- flag("admissible", NA)
    out$llReplicated <- flag("llReplicated", NA)
    bounds <- if (!is.null(diag)) diag$message[diag$type == "boundary"]
    out$boundary <- if (length(bounds) > 0) {
        paste(bounds, collapse = "; ")
    } else {
        NA_character_
    }
    out$condNumber <- flag("condNumber", NA_real_)
    out$diagnostics <- if (!is.null(diag) && nrow(diag) > 0) {
        paste(unique(diag$type), collapse = "; ")
    } else {
        NA_character_
    }

    ## Correlations between the components of the two variables
    for (f in c("trait.cor", "ar.cor", "state.cor")) {
        v <- if (!is.null(pcResults)) pcResults[[f]]
        out[[f]] <- if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v[1])
    }
    out
}
