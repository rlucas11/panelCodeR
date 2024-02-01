## General function to add constraint to parameter table
.buildConstraint <- function(parTable, lhs, op, rhs, user = 1) {
    constraint <- as.data.frame(list(
        lhs,
        op,
        rhs,
        user,
        0, 0, 0, NA, 0, "", ""
    ))
    names(constraint) <- names(parTable)
    parTable <- rbind(parTable, constraint)
    return(parTable)
}


## Constrain stability for AR part of model
.constrainStability <- function(parTable, info) {
    waves <- info$gen$maxWaves
    for (w in 3:(waves)) {
        cVals <- list(
            parTable,
            paste0("a", w),
            "==",
            paste0("a", 2)
        )
        parTable <- do.call(.buildConstraint, cVals)
    }
    if (info$gen$yVar == TRUE) {
        for (w in 3:(waves)) {
            cVals <- list(
                parTable,
                paste0("b", w),
                "==",
                paste0("b", 2)
            )
            parTable <- do.call(.buildConstraint, cVals)
        }
    }
    return(parTable)
}

## Constrain cross-lagged paths
.constrainCl <- function(parTable, info) {
    waves <- info$gen$maxWaves
    for (w in 3:waves) {
        cVals <- list(
            parTable,
            paste0("c", w),
            "==",
            paste0("c", 2)
        )
        parTable <- do.call(.buildConstraint, cVals)
    }
    for (w in 3:waves) {
        cVals <- list(
            parTable,
            paste0("d", w),
            "==",
            paste0("d", 2)
        )
        parTable <- do.call(.buildConstraint, cVals)
    }
    return(parTable)
}

## Impose stationarity on AR part of model.
.arStationarity <- function(parTable, info, constrainCor=TRUE, zero=FALSE) {
    ## Constraint for X Variance
    ## Wave 2
    if (info$gen$yVar==TRUE) {
        xvar_c <- "xvar1 - a2*a2*xvar1 - d2*d2*yvar1 - 2*a2*cov_ar1*d2"
    } else {
        xvar_c <- "xvar1 - a2*a2*xvar1"
    }
    parTable <- do.call(
        .buildConstraint,
        list(
            parTable,
            "xvar2",
            "==",
            xvar_c
        )
    )
    ## Waves beyond 2
    if (info$gen$maxWave > 2) {
        for (i in 3:info$gen$maxWave) {
            parTable <- do.call(
                .buildConstraint,
                list(
                    parTable,
                    paste0("xvar", i),
                    "==",
                    "xvar2"
                )
            )
        }
    }
    if (info$gen$yVar == TRUE) {
        ## Constraint for Y Variance
        parTable <- do.call(
            .buildConstraint,
            list(
                parTable,
                "yvar2",
                "==",
                "yvar1 - b2*b2*yvar1 - c2*c2*xvar1 - 2*b2*cov_ar1*c2"
            )
        )
        if (info$gen$maxWave > 2) {
            for (i in 3:info$gen$maxWave) {
                parTable <- do.call(
                    .buildConstraint,
                    list(
                        parTable,
                        paste0("yvar", i),
                        "==",
                        "yvar2"
                    )
                )
            }
        }
        ## Wave 2 Covariance
        if (constrainCor == TRUE) {
            if (zero == FALSE) {
                parTable <- do.call(
                    .buildConstraint,
                    list(
                        parTable,
                        "cov_ar2",
                        "==",
                        "(1-a2*b2-c2*d2)*cov_ar1-a2*c2*xvar1-b2*d2*yvar1"
                    )
                )
            } else {
                parTable <- do.call(
                    .buildConstraint,
                    list(
                        parTable,
                        "cov_ar2",
                        "==",
                        0
                    )
                )
            }
            
            ## Covariance in Waves Beyond 2
            if (info$gen$maxWave > 2) {
                for (i in 3:info$gen$maxWave) {
                    parTable <- do.call(
                        .buildConstraint,
                        list(
                            parTable,
                            paste0("cov_ar", i),
                            "==",
                            "cov_ar2"
                        )
                    )
                }
            }
        }
    }
    return(parTable)
}

## Constrain state variances
.constrainStateVar <- function(parTable, info) {
    waves <- info$gen$maxWaves
    for (w in 2:waves) {
        cVals <- list(
            parTable,
            paste0("sx", w),
            "==",
            "sx1"
        )
        parTable <- do.call(.buildConstraint, cVals)
    }
    if (info$gen$yVar == TRUE) {
        for (w in 2:waves) {
            cVals <- list(
                parTable,
                paste0("sy", w),
                "==",
                "sy1"
            )
            parTable <- do.call(.buildConstraint, cVals)
        }
    }
    ## if (xState == FALSE) {
    ##     .buildConstraint(
    ##         parTable,
    ##         sVarParams[which(sVarParams$lhs_v == info$x$name &
    ##             sVarParams$lhs_w == 1), "plabel"],
    ##         "==",
    ##         0
    ##     )
    ## }
    ## if (yState == FALSE & info$gen$yVar == TRUE) {
    ##     .buildConstraint(
    ##         parTable,
    ##         sVarParams[which(sVarParams$lhs_v == info$y$name &
    ##             sVarParams$lhs_w == 1), "plabel"],
    ##         "==",
    ##         0
    ##     )
    ## }
    return(parTable)
}

.constrainStateCors <- function(parTable, info, zero = FALSE) {
    waves <- info$gen$maxWaves
    for (w in 2:waves) {
        cVals <- list(
            parTable,
            paste0("cov_s", w),
            "==",
            "cov_s1"
        )
        parTable <- do.call(.buildConstraint, cVals)
    }
    if (zero==TRUE) {
        cVals <- list(
            parTable,
            "cov_s1",
            "==",
            0
        )
        parTable <- do.call(.buildConstraint, cVals)
    }
    return(parTable)
}

.buildLimits <- function(parTable,
                         info,
                         ar = TRUE,
                         trait = TRUE,
                         stability = TRUE,
                         cl = TRUE,
                         state = TRUE,
                         traitCors = TRUE,
                         arCors = TRUE,
                         stateCors = TRUE,
                         residCors = TRUE) {
    waves <- info$gen$maxWaves
    if (ar == TRUE) {
        parTable <- .buildConstraint(
            parTable,
            "xvar1",
            ">",
            0
        )
        if (info$gen$yVar == TRUE) {
            parTable <- .buildConstraint(
                parTable,
                "yvar1",
                ">",
                0
            )
        }
    }
    if (trait == TRUE) {
        parTable <- .buildConstraint(
            parTable,
            "x_tVar",
            ">",
            0
        )
        if (info$gen$yVar == TRUE) {
            parTable <- .buildConstraint(
                parTable,
                "y_tVar",
                ">",
                0
            )
        }
    }
    if (state == TRUE) {
        parTable <- .buildConstraint(
            parTable,
            "sx1",
            ">",
            0
        )
        if (info$gen$yVar == TRUE) {
            parTable <- .buildConstraint(
                parTable,
                "sy1",
                ">",
                0
            )
        }
        if (info$x$indicators > 1) {
            parTable <- .buildConstraint(
                parTable,
                "x1_1v",
                ">",
                0
            )
        }
        if (info$gen$yVar==TRUE) {
            if (info$y$indicators > 1) {
                parTable <- .buildConstraint(
                    parTable,
                    "y1_1v",
                    ">",
                    0
                )
            }
        }
    }
    return(parTable)
}

