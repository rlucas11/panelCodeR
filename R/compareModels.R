#' Run and Compare Models Estimated Using panelCodeR
#'
#' `compareModels()` produces code for, runs, and compares various models for
#' analyzing panel data. It produces a table of basic results. 
#'
#' @param data Dataframe with multiwave data. 
#' @param program Name of the program. Defaults to mplus as it is much quicker.
#' @param models List of models to compare. Defaults to list("starts", "riclpm",
#'   "arts", "sts", "clpm"). Can also include "dpm_c", "dpm_p", "gclm", "alt",
#'   "lgm_sm", or "lgcm" (see panelCodeR documentation). 
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
    out
}
