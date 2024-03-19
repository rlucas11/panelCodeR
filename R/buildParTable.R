.buildObserved <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar

    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }

    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character()
    )

    if (info$x$indicators > 1) {
        indList <- list(
            ustart = c(1, rep(NA, (info$x$indicators - 1))),
            free = c(0, rep(1, (info$x$indicators - 1)))
        )

        for (w in info$x$actualWaves) {
            for (i in 1:info$x$indicators) {
                loadingParTable <- data.frame(
                    lhs = paste("l", xName, w, sep = "_"),
                    op = "=~",
                    rhs = paste(xName, w, i, sep="_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = indList$free[i],
                    ustart = indList$ustart[i],
                    exo = 0,
                    label = paste0("x", w, letters[i])
                )
                initialParTable <- rbind(
                    initialParTable,
                    loadingParTable
                )
            }
        }
    } else {
        for (w in info$x$actualWaves) {
            loadingParTable <- data.frame(
                lhs = paste("l", xName, w, sep = "_"),
                op = "=~",
                rhs = paste(xName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 1,
                exo = 0,
                label = ""
            )
            var1ParTable <- data.frame(
                lhs = paste(xName, w, sep = "_"),
                op = "~~",
                rhs = paste(xName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 0,
                exo = 0,
                label = ""
            )
            initialParTable <- rbind(
                initialParTable,
                loadingParTable,
                var1ParTable
            )
        }
    }
    for (w in info$x$waves) {
        var2ParTable <- data.frame(
            lhs = paste("l", xName, w, sep = "_"),
            op = "~~",
            rhs = paste("l", xName, w, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 0,
            ustart = 0,
            exo = 0,
            label = ""
        )
        initialParTable <- rbind(
            initialParTable,
            var2ParTable
        )
    }


    if (yVar == TRUE) {
        if (info$y$indicators > 1) {
            indList <- list(
                ustart = c(1, rep(NA, (info$y$indicators - 1))),
                free = c(0, rep(1, (info$y$indicators - 1)))
            )

            for (w in info$y$actualWaves) {
                for (i in 1:info$y$indicators) {
                    loadingParTable <- data.frame(
                        lhs = paste("l", yName, w, sep = "_"),
                        op = "=~",
                        rhs = paste(yName, w, i, sep = "_"),
                        user = 1,
                        block = 1,
                        group = 1,
                        free = indList$free[i],
                        ustart = indList$ustart[i],
                        exo = 0,
                        label = paste0("y", w, letters[i])
                    )
                    initialParTable <- rbind(
                        initialParTable,
                        loadingParTable
                    )
                }
            }
        } else {
            for (w in info$y$actualWaves) {
                loadingParTable <- data.frame(
                    lhs = paste("l", yName, w, sep = "_"),
                    op = "=~",
                    rhs = paste(yName, w, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 0,
                    ustart = 1,
                    exo = 0,
                    label = ""
                )
                var1ParTable <- data.frame(
                    lhs = paste(yName, w, sep = "_"),
                    op = "~~",
                    rhs = paste(yName, w, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 0,
                    ustart = 0,
                    exo = 0,
                    label = ""
                )
                initialParTable <- rbind(
                    initialParTable,
                    loadingParTable,
                    var1ParTable
                )
            }
        }
        for (w in info$y$waves) {
            var2ParTable <- data.frame(
                lhs = paste("l", yName, w, sep = "_"),
                op = "~~",
                rhs = paste("l", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 0,
                exo = 0,
                label = ""
            )
            initialParTable <- rbind(
                initialParTable,
                var2ParTable
            )
        }

    }
    finalParTable <- initialParTable
    finalParTable$from <- "observed"
    return(finalParTable)
}

.buildPhantom <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar

    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }

    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character(),
        from = character()
    )

    phantomX <- setdiff(info$x$waves, info$x$actualWaves)
    for (i in phantomX) {
        pParTable <- data.frame(
            lhs = paste("l", xName, i, sep = "_"),
            op = "=~",
            rhs = paste("l", xName, i, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 0,
            ustart = 0,
            exo = 0,
            label = "",
            from = "phantom"
        )
        initialParTable <- rbind(
            initialParTable,
            pParTable
        )
    }
    if (yVar == TRUE) {
        phantomY <- setdiff(info$y$waves, info$y$actualWaves)
        for (i in phantomY) {
            pParTable <- data.frame(
                lhs = paste("l", yName, i, sep = "_"),
                op = "=~",
                rhs = paste("l", yName, i, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 0,
                exo = 0,
                label = "",
                from = "phantom"
            )
            initialParTable <- rbind(
                initialParTable,
                pParTable
            )
        }
    }
    return(initialParTable)
}


    
.buildResidVar <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    
    xName <- info$x$name
    xInd <- info$x$indicators
    xWaves <- info$x$actualWaves
    if (yVar == TRUE) {
        yName <- info$y$name
        yInd <- info$y$indicators
        yWaves <- info$y$actualWaves
    }

    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character(),
        from = character()
    )

    if (xInd > 1) {
        for (i in 1:xInd) {
            for (j in xWaves) {
                rCorTable <- data.frame(
                    lhs = paste(xName, j, i, sep = "_"),
                    op = "~~",
                    rhs = paste(xName, j, i, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = paste0("x",
                                   j,
                                   "_",
                                   i,
                                   "v"),
                    from = "residVar"
                )
                initialParTable <- rbind(
                    initialParTable,
                    rCorTable
                )
            }
        }
    }


    if (yVar == TRUE) {
        if(yInd > 1) {
            for (i in 1:yInd) {
                for (j in yWaves) {
                    rCorTable <- data.frame(
                        lhs = paste(yName, j, i, sep = "_"),
                        op = "~~",
                        rhs = paste(yName, j, i, sep = "_"),
                        user = 1,
                        block = 1,
                        group = 1,
                        free = 1,
                        ustart = NA,
                        exo = 0,
                        label = paste0("y",
                                       j,
                                       "_",
                                       i,
                                       "v"),
                        from = "residVar"
                    )
                    initialParTable <- rbind(
                        initialParTable,
                        rCorTable
                    )
                }
            }
        }
    }
    finalParTable <- initialParTable[order(initialParTable$op), ]
    return(finalParTable)  
}


.buildAr <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    
    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }

    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character()
    ) 
    
    for (w in info$x$waves) {
        loadingParTable <- data.frame(
            lhs = paste("a", xName, w, sep="_"),
            op = "=~",
            rhs = paste("l", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 0,
            ustart = 1,
            exo = 0,
            label = ""
        )
        var1ParTable <- data.frame(
            lhs = paste("a", xName, w, sep = "_"),
            op = "~~",
            rhs = paste("a", xName, w, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = paste0("xvar", w)
        )
        initialParTable <- rbind(
            initialParTable,
            loadingParTable,
            var1ParTable
        )
    }

    if (yVar == TRUE) {
        for (w in info$y$waves) {
            loadingParTable <- data.frame(
                lhs = paste("a", yName, w, sep="_"),
                op = "=~",
                rhs = paste("l", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 1,
                exo = 0,
                label = ""
            )
            var1ParTable <- data.frame(
                lhs = paste("a", yName, w, sep = "_"),
                op = "~~",
                rhs = paste("a", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("yvar", w)
            )
            initialParTable <- rbind(
                initialParTable,
                loadingParTable,
                var1ParTable
            )
        }
    }
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "ar"
    return(finalParTable)  
}



.buildTrait <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar

    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }


    
    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character()
    )
    
    for (w in info$x$waves) {
        loadingParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "=~",
            rhs = paste("l", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 0,
            ustart = 1,
            exo = 0,
            label = ""
        )
        initialParTable <- rbind(
            initialParTable,
            loadingParTable
        )
    }

    var1ParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "~~",
            rhs = paste("t", xName, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "x_tVar"
    )
    initialParTable <- rbind(
        initialParTable,
        var1ParTable
    )

    if (yVar == TRUE) {
        for (w in info$y$waves) {
            loadingParTable <- data.frame(
                lhs = paste("t", yName, sep = "_"),
                op = "=~",
                rhs = paste("l", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 1,
                exo = 0,
                label = ""
            )
            initialParTable <- rbind(
                initialParTable,
                loadingParTable
            )
        }

        var1ParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("t", yName, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "y_tVar"
        )
        initialParTable <- rbind(
            initialParTable,
            var1ParTable
        )
    }
    
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "trait"
    return(finalParTable)  
}


.buildDpmTrait <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar

    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }

    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character()
    )
    
    for (w in info$x$waves[-1]) {
        loadingParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "=~",
            rhs = paste("a", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 0,
            ustart = 1,
            exo = 0,
            label = ""
        )
        initialParTable <- rbind(
            initialParTable,
            loadingParTable
        )
    }

    var1ParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "~~",
            rhs = paste("t", xName, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "x_tVar"
    )
    initialParTable <- rbind(
        initialParTable,
        var1ParTable
    )

    if (yVar == TRUE) {
        for (w in info$y$waves[-1]) {
            loadingParTable <- data.frame(
                lhs = paste("t", yName, sep = "_"),
                op = "=~",
                rhs = paste("a", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 1,
                exo = 0,
                label = ""
            )
            initialParTable <- rbind(
                initialParTable,
                loadingParTable
            )
        }

        var1ParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("t", yName, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "y_tVar"
        )
        initialParTable <- rbind(
            initialParTable,
            var1ParTable
        )
    }
    
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "trait"
    return(finalParTable)  
}



.buildStability <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar

    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }


    
    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character()
    )
    
    for (w in info$x$waves[-1]) {
        stabParTable <- data.frame(
            lhs = paste("a", xName, w, sep = "_"),
            op = "~",
            rhs = paste("a", xName, (w - 1), sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = paste0("a", w)
        )
        initialParTable <- rbind(
            initialParTable,
            stabParTable
        )
    }

    for (w in info$y$waves[-1]) {
        stabParTable <- data.frame(
            lhs = paste("a", yName, w, sep = "_"),
            op = "~",
            rhs = paste("a", yName, (w - 1), sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = paste0("b", w)
        )
        initialParTable <- rbind(
            initialParTable,
            stabParTable
        )
    }
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "stability"
    return(finalParTable)  
}

.buildCrossLag <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    
    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }

    if (yVar == TRUE) {
        ## Create initial table
        initialParTable <- data.frame(
            lhs = character(),
            op = character(),
            rhs = character(),
            user = integer(),
            block = integer(),
            group = integer(),
            free = integer(),
            ustart = numeric(),
            exo = integer(),
            label = character()
        )

        for (w in info$x$waves[-1]) {
            cl1ParTable <- data.frame(
                lhs = paste("a", xName, w, sep = "_"),
                op = "~",
                rhs = paste("a", yName, (w - 1), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("d", w)
            )
            cl2ParTable <- data.frame(
                lhs = paste("a", yName, w, sep = "_"),
                op = "~",
                rhs = paste("a", xName, (w - 1), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("c", w)
            )
            initialParTable <- rbind(
                initialParTable,
                cl1ParTable,
                cl2ParTable
            )
        }
        finalParTable <- initialParTable[order(initialParTable$op), ]
        finalParTable$from <- "cl"
        return(finalParTable)
    }
}

.buildState <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    
    xName <- info$x$name
    if (yVar == TRUE) {
        yName <- info$y$name
    }


    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character()
    )
    
    for (w in info$x$waves) {
        loadingParTable <- data.frame(
            lhs = paste("s", xName, w, sep="_"),
            op = "=~",
            rhs = paste("l", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 0,
            ustart = 1,
            exo = 0,
            label = ""
        )
        var1ParTable <- data.frame(
            lhs = paste("s", xName, w, sep = "_"),
            op = "~~",
            rhs = paste("s", xName, w, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = paste0("sx", w)
        )
        initialParTable <- rbind(
            initialParTable,
            loadingParTable,
            var1ParTable
        )
    }

    if (yVar == TRUE) {
        for (w in info$y$waves) {
            loadingParTable <- data.frame(
                lhs = paste("s", yName, w, sep="_"),
                op = "=~",
                rhs = paste("l", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 0,
                ustart = 1,
                exo = 0,
                label = ""
            )
            var1ParTable <- data.frame(
                lhs = paste("s", yName, w, sep = "_"),
                op = "~~",
                rhs = paste("s", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("sy", w)
            )
            initialParTable <- rbind(
                initialParTable,
                loadingParTable,
                var1ParTable
            )
        }
    }
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "state"
    return(finalParTable)  
}


.buildCors <- function(info, ar = TRUE, trait = TRUE, state = TRUE) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    xName <- info$x$name
    
    if (yVar == TRUE) {
        yName <- info$y$name

        if (trait==TRUE) {
            corParTable <- data.frame(
                lhs = paste("t", xName, sep = "_"),
                op = "~~",
                rhs = paste("t", yName, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "cov_txty"
            ) 
        } else {
            corParTable <- NULL
        }
        
        if (ar == TRUE) {
            ## Create initial table
            arCorParTable <- data.frame(
                lhs = character(),
                op = character(),
                rhs = character(),
                user = integer(),
                block = integer(),
                group = integer(),
                free = integer(),
                ustart = numeric(),
                exo = integer(),
                label = character()
            )
            for (w in 1:info$gen$maxWaves) {
                newCorParTable <- data.frame(
                    lhs = paste("a", xName, w, sep = "_"),
                    op = "~~",
                    rhs = paste("a", yName, w, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = paste0("cov_ar", w)
                )
                arCorParTable <- rbind(
                    arCorParTable,
                    newCorParTable
                )
            }
        } else {
            arCorParTable <- NULL
        }
        
        if (state == TRUE) {
            stateCorParTable <- data.frame(
                lhs = character(),
                op = character(),
                rhs = character(),
                user = integer(),
                block = integer(),
                group = integer(),
                free = integer(),
                ustart = numeric(),
                exo = integer(),
                label = character()
            )
            for (w in 1:info$gen$maxWaves) {
                newStateCorParTable <- data.frame(
                    lhs = paste("s", xName, w, sep = "_"),
                    op = "~~",
                    rhs = paste("s", yName, w, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = paste0("cov_s", w)
                )
                stateCorParTable <- rbind(
                    stateCorParTable,
                    newStateCorParTable
                )
            }
        } else {
            stateCorParTable <- NULL
        }
        corParTable <- do.call(
            rbind,
            list(
                corParTable,
                arCorParTable,
                stateCorParTable
            )
        )
        corParTable <- corParTable[order(corParTable$lhs), ]
        corParTable$from <- "cors"
        return(corParTable)
    }
}

.buildDpmCors <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    xName <- info$x$name

    corParTable <- data.frame(
        lhs = paste("t", xName, sep = "_"),
        op = "~~",
        rhs = paste("a", xName, 1, sep = "_"),
        user = 1,
        block = 1,
        group = 1,
        free = 1,
        ustart = NA,
        exo = 0,
        label = "cov_TArX"
    )
    
    if (yVar == TRUE) {
        yName <- info$y$name
        ytCorParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("a", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "cov_TArY"
        )
        xyCorParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "~~",
            rhs = paste("a", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "cov_TxAy"
        )
        yxCorParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("a", xName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "cov_TyAx"
        )
        corParTable <- do.call(
        rbind,
        list(
            corParTable,
            ytCorParTable,
            xyCorParTable,
            yxCorParTable
        )
    )
    }
    
    corParTable$from <- "dpmCors"
    return(corParTable)
}


.buildResidCors <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    
    xName <- info$x$name
    xInd <- info$x$indicators
    xWaves <- info$x$actualWaves
    if (yVar == TRUE) {
        yName <- info$y$name
        yInd <- info$y$indicators
        yWaves <- info$y$actualWaves
    }

    ## Create initial table
    initialParTable <- data.frame(
        lhs = character(),
        op = character(),
        rhs = character(),
        user = integer(),
        block = integer(),
        group = integer(),
        free = integer(),
        ustart = numeric(),
        exo = integer(),
        label = character(),
        from = character()
    )

    if (xInd > 1) {
        for (i in 1:xInd) {
            for (j in 1:(length(xWaves)-1)) {
                jWave <- xWaves[j]
                for (k in xWaves[(j + 1):length(xWaves)]) {
                    rCorTable <- data.frame(
                        lhs = paste(xName, jWave, i, sep = "_"),
                        op = "~~",
                        rhs = paste(xName, k, i, sep = "_"),
                        user = 1,
                        block = 1,
                        group = 1,
                        free = 1,
                        ustart = NA,
                        exo = 0,
                        label = paste0("x",
                                       k,
                                       letters[i],
                                       "_l",
                                       k-jWave),
                        from = "residCors"
                    )
                    initialParTable <- rbind(
                        initialParTable,
                        rCorTable
                    )
                }
            }
        }
    }

    if (yVar == TRUE) {
        if(yInd > 1) {
            for (i in 1:yInd) {
                for (j in 1:(length(yWaves)-1)) {
                    jWave <- yWaves[j]
                    for (k in yWaves[(j + 1):length(yWaves)]) {
                        rCorTable <- data.frame(
                            lhs = paste(yName, jWave, i, sep = "_"),
                            op = "~~",
                            rhs = paste(yName, k, i, sep = "_"),
                            user = 1,
                            block = 1,
                            group = 1,
                            free = 1,
                            ustart = NA,
                            exo = 0,
                            label = paste0("y",
                                           k,
                                           letters[i],
                                           "_l",
                                           k-jWave),
                            from = "residCors"
                        )
                        initialParTable <- rbind(
                            initialParTable,
                            rCorTable
                        )
                    }
                }
            }
        }
    }
    return(initialParTable)
}

.buildTable <- function(info,
                        ar = TRUE,
                        trait = TRUE,
                        state = TRUE,
                        stability = TRUE,
                        cl = TRUE,
                        traitCors = TRUE,
                        arCors = TRUE,
                        stateCors = TRUE,
                        residCors = TRUE,
                        dpm = FALSE) {

    if (dpm == TRUE) {
        trait <- FALSE
        state <- FALSE
    }

    if (trait == FALSE) {
        traitCors == FALSE
    }

    if (ar == FALSE) {
        arCors <- FALSE
    }

    if (state == FALSE) {
        stateCors <- FALSE
    }
    
    components <- list(
        ph = .buildPhantom(info),
        obs = .buildObserved(info),
        ar = .buildAr(info),
        residVar = .buildResidVar(info),
        trait = .buildTrait(info),
        dpmTrait = .buildDpmTrait(info),
        stability = .buildStability(info),
        cl = .buildCrossLag(info),
        state = .buildState(info),
        cors = .buildCors(info, ar = arCors, trait = traitCors, state = stateCors),
        dpmCors = .buildDpmCors(info),
        residCors = .buildResidCors(info)
    )
    
    if (ar == FALSE) {
        components$ar <- NULL
    }
    ## Fix resid var stuff here
    if (trait == FALSE) {
        components$trait <- NULL
    }
    if (dpm == FALSE) {
        components$dpmTrait <- NULL
    }
    if (stability == FALSE) {
        components$stability <- NULL
    }
    if (cl == FALSE) {
        components$cl <- NULL
    }
    if (state == FALSE) {
        components$state <- NULL
    }
    if (residCors == FALSE) {
        components$residCors <- NULL
    }
    if (dpm == FALSE) {
        components$dpmCors <- NULL
    }
           
    
    model <- do.call(
        rbind,
        components
    )

    ## Get Unique names
    allVarNames <- unique(model$lhs)
    varNamesInfo <- data.frame(
        varName = allVarNames,
        varName1 = stringr::str_split_i(allVarNames, "_", 1),
        varName2 = stringr::str_split_i(allVarNames, "_", 2),
        varName3 = stringr::str_split_i(allVarNames, "_", 3)
    )

    ## Get component, variable, and wave info for LHS and RHS
    latentVarInfo <- varNamesInfo[which(varNamesInfo$varName1 == "l" |
        varNamesInfo$varName1 == "a" |
        varNamesInfo$varName1 == "t" |
        varNamesInfo$varName1 == "s"), ]
    names(latentVarInfo) <- c("varName", "component", "variable", "wave")
    latentVarInfo$indicator <- NA
    if(info$gen$yVar==TRUE) {
        observedVarInfo <- varNamesInfo[which(varNamesInfo$varName1 == info$x$name |
                                              varNamesInfo$varName1 == info$y$name), ]
    } else {
        observedVarInfo <- varNamesInfo[which(varNamesInfo$varName1 == info$x$name),]
    }
    names(observedVarInfo) <- c("varName", "variable", "wave", "indicator")
    observedVarInfo$component <- NA
    varInfo <- rbind(
        latentVarInfo,
        observedVarInfo[, c(
            "varName",
            "component",
            "variable",
            "wave",
            "indicator"
        )]
    )
    model <- dplyr::left_join(model, varInfo, by = c("lhs" = "varName"))
    names(model)[12:15] <- c("lhs_c", "lhs_v", "lhs_w", "lhs_i")
    model <- dplyr::left_join(model, varInfo, by = c("rhs" = "varName"))
    names(model)[16:19] <- c("rhs_c", "rhs_v", "rhs_w", "rhs_i")

    ## Create explicit intercepts
    int.ov <- data.frame(
        lhs = observedVarInfo$varName,
        op = rep("~1", length(observedVarInfo$varName)),
        rhs = rep("", length(observedVarInfo$varName)),
        user = rep(1, length(observedVarInfo$varName)),
        block = rep(1, length(observedVarInfo$varName)),
        group = rep(1, length(observedVarInfo$varName)),
        free = rep(1, length(observedVarInfo$varName)),
        ustart = rep(NA, length(observedVarInfo$varName)),
        exo = rep(0, length(observedVarInfo$varName)),
        label = rep("", length(observedVarInfo$varName)),
        from = rep("int.ov", length(observedVarInfo$varName))
    )

    int.lv <- data.frame(
        lhs = latentVarInfo$varName,
        op = rep("~1", length(latentVarInfo$varName)),
        rhs = rep("", length(latentVarInfo$varName)),
        user = rep(1, length(latentVarInfo$varName)),
        block = rep(1, length(latentVarInfo$varName)),
        group = rep(1, length(latentVarInfo$varName)),
        free = rep(0, length(latentVarInfo$varName)),
        ustart = rep(0, length(latentVarInfo$varName)),
        exo = rep(0, length(latentVarInfo$varName)),
        label = rep("", length(latentVarInfo$varName)),
        from = rep("int.lv", length(latentVarInfo$varName))
        )

    keepVars <- c(
        "lhs", "op", "rhs", "free", "label",
        paste0(
            rep(c("lhs", "rhs"), each = 4),
            c("_c", "_v", "_w", "_i")
        )
    )
              
    ## Get Free parameters
    freeParams <- model[which(model$free > 0), keepVars] 

    ## Get Variance Parameters
    varParams <- subset(freeParams, lhs == rhs)

    ## Get Stability Parameters
    stabParams <- subset(freeParams,
                         op == "~" &
                         lhs_v == rhs_v &
                         lhs_c == "a"
                         )

    ## Get Cross-Lag Parameters
    clParams <- subset(
        freeParams,
        op == "~" &
        lhs_v != rhs_v &
        lhs_c == "a"
    ) %>%
        dplyr::arrange(lhs_v, as.numeric(lhs_w))
    
    ## Get Correlation Parameters
    corParams <- subset(freeParams,
                        op == "~~" &
                        lhs_v != rhs_v
                        ) %>%
        dplyr::arrange(rhs_c, lhs_v, as.numeric(lhs_w))

    ## Get Residual Var/Cov Parameters
    residParams <- subset(
        freeParams,
        op == "~~" &
        !is.na(lhs_i)
    )

    
    ## Revert model back to original
    model <- model[c(
        "lhs",
        "op",
        "rhs",
        "user",
        "block",
        "group",
        "free",
        "ustart",
        "exo",
        "label"
    )]

    model <- rbind(model, int.ov[,1:10], int.lv[,1:10])
    ## Create plabels
    model$plabel <- paste0(".p", 1:nrow(model), ".")

    ## Combine Info
    tableInfo <- list(
        model = model,
        freeParams = freeParams,
        varParams = varParams,
        stabParams = stabParams,
        clParams = clParams,
        corParams = corParams,
        residParams = residParams
    )
    return(tableInfo)
}
