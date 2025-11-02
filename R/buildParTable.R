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
        ## Constrain to 0 after adding impulses
        var1ParTable <- data.frame(
            lhs = paste("a", xName, w, sep = "_"),
            op = "~~",
            rhs = paste("a", xName, w, sep = "_"),
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
            ## Constrain to zero after adding impulses
            var1ParTable <- data.frame(
                lhs = paste("a", yName, w, sep = "_"),
                op = "~~",
                rhs = paste("a", yName, w, sep = "_"),
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
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "ar"
    return(finalParTable)  
}

.buildImpulses <- function(info) {
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
            lhs = paste("i", xName, w, sep="_"),
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
        var1ParTable <- data.frame(
            lhs = paste("i", xName, w, sep = "_"),
            op = "~~",
            rhs = paste("i", xName, w, sep = "_"),
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
                lhs = paste("i", yName, w, sep="_"),
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
            var1ParTable <- data.frame(
                lhs = paste("i", yName, w, sep = "_"),
                op = "~~",
                rhs = paste("i", yName, w, sep = "_"),
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
    finalParTable$from <- "imp"
    return(finalParTable)  
}



.buildTrait <- function(info, free=FALSE) {
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

    ## Create loadings
    if (free == FALSE) {
        freeValues <- rep(0, length(info$x$waves))
        ustartValues <- rep(1, length(info$x$waves))
        xLabelsValues <- rep("", length(info$x$waves))
        if (yVar == TRUE) {
            yLabelsValues <- rep("", length(info$y$waves))
        }
    } else {
        freeValues <- c(0, rep(1, (length(info$x$waves)-1)))
        ustartValues <- c(1, rep(NA, (length(info$x$waves)-1)))
        xLabelsValues <- paste0("xtl_", info$x$waves)
        if (yVar == TRUE) {
            yLabelsValues <- paste0("ytl_", info$y$waves)
        }
    }

    
    
    
    for (w in info$x$waves) {
        loadingParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "=~",
            rhs = paste("l", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = freeValues[w],
            ustart = ustartValues[w],
            exo = 0,
            label = xLabelsValues[w]
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
                free = freeValues[w],
                ustart = ustartValues[w],
                exo = 0,
                label = yLabelsValues[w]
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

#' @importFrom stats median
.buildGclmTrait <- function(info, free=TRUE) {
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

    ## Create loadings
    if (free == FALSE) {
        freeValues <- rep(0, length(info$x$waves))
        ustartValues <- rep(1, length(info$x$waves))
        xLabelsValues <- rep("", length(info$x$waves))
        if (yVar == TRUE) {
            yLabelsValues <- rep("", length(info$y$waves))
        }
    } else {
        freeValues <- c(0, rep(1, (length(info$x$waves)-1)))
        ustartValues <- c(1, rep(NA, (length(info$x$waves)-1)))
        xLabelsValues <- paste0("xtl_", info$x$waves)
        if (yVar == TRUE) {
            yLabelsValues <- paste0("ytl_", info$y$waves)
        }
    }
    
    for (w in info$x$waves) {
        loadingParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "=~",
            rhs = paste("a", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = freeValues[w],
            ustart = ustartValues[w],
            exo = 0,
            label = xLabelsValues[w]
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
                rhs = paste("a", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = freeValues[w],
                ustart = ustartValues[w],
                exo = 0,
                label = yLabelsValues[w]
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


.buildSlope <- function(info, free=FALSE, slope="none") {
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

    ## Create loadings
    if (slope == "linear") {
        freeValues <- rep(0, length(info$x$waves))
        ustartValues <- seq(from = 0, to = length(info$x$waves))
        xLabelsValues <- rep("", length(info$x$waves))
        if (yVar == TRUE) {
            yLabelsValues <- rep("", length(info$y$waves))
        }
    } else if (slope == "centered") {
        firstLoading <- (1 - median(1:info$gen$maxWaves))
        freeValues <- rep(0, length(info$x$waves))
        ustartValues <- seq(from = firstLoading,
                            to = firstLoading + info$gen$maxWaves)
        xLabelsValues <- rep("", length(info$x$waves))
        if (yVar == TRUE) {
            yLabelsValues <- rep("", length(info$y$waves))
        }
    } else if (slope == "basis") {
        freeValues <- c(0, rep(1, (length(info$x$waves)-2)), 0)
        ustartValues <- c(0, rep(NA, (length(info$x$waves)-2)), 1)
        xLabelsValues <- rep("", length(info$x$waves))
        if (yVar == TRUE) {
            yLabelsValues <- rep("", length(info$y$waves))
        }
    } else
        return(NULL)

    for (w in info$x$waves) {
        loadingParTable <- data.frame(
            lhs = paste("sl", xName, sep = "_"),
            op = "=~",
            rhs = paste("l", xName, w, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = freeValues[w],
            ustart = ustartValues[w],
            exo = 0,
            label = xLabelsValues[w]
        )
        initialParTable <- rbind(
            initialParTable,
            loadingParTable
        )
    }

    var1ParTable <- data.frame(
            lhs = paste("sl", xName, sep = "_"),
            op = "~~",
            rhs = paste("sl", xName, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "x_slVar"
    )
    initialParTable <- rbind(
        initialParTable,
        var1ParTable
    )

    if (yVar == TRUE) {
        for (w in info$y$waves) {
            loadingParTable <- data.frame(
                lhs = paste("sl", yName, sep = "_"),
                op = "=~",
                rhs = paste("l", yName, w, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = freeValues[w],
                ustart = ustartValues[w],
                exo = 0,
                label = yLabelsValues[w]
            )
            initialParTable <- rbind(
                initialParTable,
                loadingParTable
            )
        }

        var1ParTable <- data.frame(
            lhs = paste("sl", yName, sep = "_"),
            op = "~~",
            rhs = paste("sl", yName, sep="_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "y_slVar"
        )
        initialParTable <- rbind(
            initialParTable,
            var1ParTable
        )
    }
    
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "slope"
    return(finalParTable)  
}


.buildStability <- function(info, lags) {
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

    for (l in 1:lags) {
        for (w in info$x$waves[-c(1:l)]) {
            pLabel <- paste("a", l, w, sep = "_")
            stabParTable <- data.frame(
                lhs = paste("a", xName, w, sep = "_"),
                op = "~",
                rhs = paste("a", xName, (w - l), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = pLabel
            )
            initialParTable <- rbind(
                initialParTable,
                stabParTable
            )
        }
    }

    for (l in 1:lags) {
        for (w in info$y$waves[-c(1:l)]) {
            pLabel <- paste("b", l, w, sep = "_")
            stabParTable <- data.frame(
                lhs = paste("a", yName, w, sep = "_"),
                op = "~",
                rhs = paste("a", yName, (w - l), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = pLabel
            )
            initialParTable <- rbind(
                initialParTable,
                stabParTable
            )
        }
    }
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "stability"
    return(finalParTable)  
}




.buildCrossLag <- function(info, lags) {
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

        for (l in 1:lags) {
            for (w in info$x$waves[-c(1:l)]) {
                xpLabel <- paste("d", l, w, sep = "_")
                ypLabel <- paste("c", l, w, sep = "_")
                cl1ParTable <- data.frame(
                    lhs = paste("a", xName, w, sep = "_"),
                    op = "~",
                    rhs = paste("a", yName, (w - l), sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = xpLabel
                )
                cl2ParTable <- data.frame(
                    lhs = paste("a", yName, w, sep = "_"),
                    op = "~",
                    rhs = paste("a", xName, (w - l), sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = ypLabel
                )
                initialParTable <- rbind(
                    initialParTable,
                    cl1ParTable,
                    cl2ParTable
                )
            }
        }
        finalParTable <- initialParTable[order(initialParTable$op), ]
        finalParTable$from <- "cl"
        return(finalParTable)
    }
}


.buildMa <- function(info) {
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
        maParTable <- data.frame(
            lhs = paste("a", xName, w, sep = "_"),
            op = "~",
            rhs = paste("i", xName, (w - 1), sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = paste0("ma_a", w)
        )
        initialParTable <- rbind(
            initialParTable,
            maParTable
        )
    }

    if (yVar == TRUE) {
        for (w in info$y$waves[-1]) {
            maParTable <- data.frame(
                lhs = paste("a", yName, w, sep = "_"),
                op = "~",
                rhs = paste("i", yName, (w - 1), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("ma_b", w)
            )
            initialParTable <- rbind(
                initialParTable,
                maParTable
            )
        }
    }
    finalParTable <- initialParTable[order(initialParTable$op), ]
    finalParTable$from <- "stability"
    return(finalParTable)  
}


.buildClMa <- function(info) {
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
            ma1ParTable <- data.frame(
                lhs = paste("a", xName, w, sep = "_"),
                op = "~",
                rhs = paste("i", yName, (w - 1), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("ma_d", w)
            )
            ma2ParTable <- data.frame(
                lhs = paste("a", yName, w, sep = "_"),
                op = "~",
                rhs = paste("i", xName, (w - 1), sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = paste0("ma_c", w)
            )
            initialParTable <- rbind(
                initialParTable,
                ma1ParTable,
                ma2ParTable
            )
        }
        finalParTable <- initialParTable[order(initialParTable$op), ]
        finalParTable$from <- "ma"
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


.buildCors <- function(info,
                       ar = TRUE,
                       trait = TRUE,
                       state = TRUE,
                       slope = "none") {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    xName <- info$x$name
    
    if (yVar == TRUE) {
        yName <- info$y$name

        if (trait == TRUE) {
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

        if (slope != "none") {
            slopeCorParTable <- data.frame(
                lhs = c(
                    paste("t", xName, sep = "_"),
                    paste("t", yName, sep = "_"),
                    paste("sl", xName, sep = "_")
                ),
                op = rep("~~",3),
                rhs = c(
                    paste("sl", xName, sep = "_"),
                    paste("sl", yName, sep = "_"),
                    paste("sl", yName, sep = "_")
                ),
                user = rep(1, 3),
                block = rep(1, 3),
                group = rep(1, 3),
                free = rep(1, 3),
                ustart = rep(NA, 3),
                exo = rep(0, 3),
                label = c(
                    "cov_txsx",
                    "cov_tysy",
                    "cov_sxsy"
                    )
            )
        } else {
            slopeCorParTable <- NULL
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
                    lhs = paste("i", xName, w, sep = "_"),
                    op = "~~",
                    rhs = paste("i", yName, w, sep = "_"),
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
                slopeCorParTable,
                arCorParTable,
                stateCorParTable
            )
        )
        if (!is.null(corParTable)) {
            corParTable <- corParTable[order(corParTable$lhs), ]
            corParTable$from <- "cors"
        }
        return(corParTable)
    } else {
        if (slope != "none") {
            slopeCorParTable <- data.frame(
                lhs = paste("t", xName, sep = "_"),
                op = "~~",
                rhs = paste("sl", xName, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "cov_txsx"
            )
            corParTable <- slopeCorParTable
            corParTable <- corParTable[order(corParTable$lhs), ]
            corParTable$from <- "cors"
        } else {
            corParTable <- NULL
        }
        return(corParTable)
    }
}

.buildPredetermined <- function(info, slope) {
    yVar <- info$gen$yVar
    xName <- info$x$name

    if (yVar == TRUE) {
        yName <- info$y$name
    }

    corParTable <- data.frame(
        lhs = paste("t", xName, sep = "_"),
        op = "~~",
        rhs = paste("i", xName, 1, sep = "_"),
        user = 1,
        block = 1,
        group = 1,
        free = 1,
        ustart = NA,
        exo = 0,
        label = "c_tx_ix"
    )

    if (slope != "none") {
        sCorParTable <- data.frame(
            lhs = paste("sl", xName, sep = "_"),
            op = "~~",
            rhs = paste("i", xName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_slx_ix"
        )
        corParTable <- do.call(
            rbind,
            list(
                corParTable,
                sCorParTable
            )
        )
    }

    if (yVar == TRUE) {
        yCorParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("i", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_ty_iy"
        )
        y2CorParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "~~",
            rhs = paste("i", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_tx_iy"
        )
        y3CorParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("i", xName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_ty_ix"
        )

        corParTable <- do.call(
            rbind,
            list(
                corParTable,
                yCorParTable,
                y2CorParTable,
                y3CorParTable
            )
        )

        if (slope != "none") {
            syCorParTable <- data.frame(
                lhs = paste("sl", yName, sep = "_"),
                op = "~~",
                rhs = paste("i", yName, 1, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_sly_ix"
            )
            sy2CorParTable <- data.frame(
                lhs = paste("sl", xName, sep = "_"),
                op = "~~",
                rhs = paste("i", yName, 1, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_slx_iy"
            )
            corParTable <- do.call(
                rbind,
                list(
                    corParTable,
                    syCorParTable,
                    sy2CorParTable
                )
            )
        }
    }
    corParTable$from <- "cors"
    return(corParTable)
}


.buildMeasurement <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    xName <- info$x$name
    waves <- info$gen$maxWaves
    xiNames <- paste(
        rep("i", waves),
        rep(xName, waves),
        1:waves,
        sep="_"
        )
    
    if (yVar == TRUE) {
        yName <- info$y$name
        yiNames <- paste(
            rep("i", waves),
            rep(yName, waves),
            1:waves,
            sep="_"
        )
        iNames <- c(xiNames, yiNames)
    } else {
        iNames <- xiNames
    }
    
    ## Create initial table
    corParTable <- data.frame(
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
    for (i in 1:(length(iNames)-1)) {
        for (j in (i+1):length(iNames)) {
            newCorParTable <- data.frame(
                lhs = iNames[i],
                op = "~~",
                rhs = iNames[j],
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = ""
            )
            corParTable <- rbind(
                corParTable,
                newCorParTable
            )
        }
    }
    corParTable$from <- "cors"
    return(corParTable)
} 

.buildDpmTraitCors <- function(info) {

    yVar <- info$gen$yVar
    xName <- info$x$name

    if (yVar == TRUE) {
        yName <- info$y$name
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
        corParTable$from <- "dpmTraitCors"
        return(corParTable)
    } else {
        return(NULL)
    }
}
  

.buildDpmCors <- function(info, slope) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    xName <- info$x$name

    corParTable <- data.frame(
        lhs = paste("t", xName, sep = "_"),
        op = "~~",
        rhs = paste("i", xName, 1, sep = "_"),
        user = 1,
        block = 1,
        group = 1,
        free = 1,
        ustart = NA,
        exo = 0,
        label = "c_tx_ix"
    )

    if (slope != "none") {
        ## slCorParTable <- data.frame(
        ##     lhs = paste("t", xName, sep = "_"),
        ##     op = "~~",
        ##     rhs = paste("sl", xName, sep = "_"),
        ##     user = 1,
        ##     block = 1,
        ##     group = 1,
        ##     free = 1,
        ##     ustart = NA,
        ##     exo = 0,
        ##     label = "c_tx_slx"
        ## )
        ## Removed because redundant
        slArCorParTable <- data.frame(
            lhs = paste("sl", xName, sep = "_"),
            op = "~~",
            rhs = paste("i", xName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_slx_ix"
        )
        corParTable <- do.call(
            rbind,
            list(
                corParTable,
##                slCorParTable,
                slArCorParTable
            )
        )
    }
    
    
    if (yVar == TRUE) {
        yName <- info$y$name
        y1CorParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("i", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_ty_iy"
        )
        y2CorParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "~~",
            rhs = paste("i", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_tx_iy"
        )
        y3CorParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "~~",
            rhs = paste("i", xName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "c_ty_ix"
        )
        corParTable <- do.call(
            rbind,
            list(
                corParTable,
                y1CorParTable,
                y2CorParTable,
                y3CorParTable
            )
        )
        if (slope != "none") {
            ## ysCorParTable <- data.frame(
            ##     lhs = paste("sl", yName, sep = "_"),
            ##     op = "~~",
            ##     rhs = paste("t", xName, sep = "_"),
            ##     user = 1,
            ##     block = 1,
            ##     group = 1,
            ##     free = 1,
            ##     ustart = NA,
            ##     exo = 0,
            ##     label = "c_sly_ty"
            ## )
            ## Removed because redundant
            ys2CorParTable <- data.frame(
                lhs = paste("sl", yName, sep = "_"),
                op = "~~",
                rhs = paste("i", yName, 1, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_sly_iy"
            )
            ## ys3CorParTable <- data.frame(
            ##     lhs = paste("sl", yName, sep = "_"),
            ##     op = "~~",
            ##     rhs = paste("sl", xName, sep = "_"),
            ##     user = 1,
            ##     block = 1,
            ##     group = 1,
            ##     free = 1,
            ##     ustart = NA,
            ##     exo = 0,
            ##     label = "c_sly_slx"
            ## )
            ys4CorParTable <- data.frame(
                lhs = paste("sl", yName, sep = "_"),
                op = "~~",
                rhs = paste("i", xName, 1, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_sly_ix"
            )
            ys5CorParTable <- data.frame(
                lhs = paste("t", yName, sep = "_"),
                op = "~~",
                rhs = paste("sl", xName, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_ty_slx"
            )
            ys6CorParTable <- data.frame(
                lhs = paste("t", xName, sep = "_"),
                op = "~~",
                rhs = paste("sl", yName, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_tx_sly"
            )
            ys7CorParTable <- data.frame(
                lhs = paste("sl", xName, sep = "_"),
                op = "~~",
                rhs = paste("i", yName, 1, sep = "_"),
                user = 1,
                block = 1,
                group = 1,
                free = 1,
                ustart = NA,
                exo = 0,
                label = "c_slx_iy"
            )
            corParTable <- do.call(
                rbind,
                list(
                    corParTable,
                ##    ysCorParTable,
                    ys2CorParTable,
                ##    ys3CorParTable,
                    ys4CorParTable,
                    ys5CorParTable,
                    ys6CorParTable,
                    ys7CorParTable
                )
            )
        }
    }
    
    corParTable$from <- "dpmCors"
    return(corParTable)
}

.buildDpmLoadings <- function(info) {
    ## Check if bivariate or univariate
    yVar <- info$gen$yVar
    xName <- info$x$name

    loadParTable <- data.frame(
        lhs = paste("t", xName, sep = "_"),
        op = "=~",
        rhs = paste("a", xName, 1, sep = "_"),
       user = 1,
        block = 1,
        group = 1,
        free = 1,
        ustart = NA,
        exo = 0,
        label = "d_a"
    )
    
    if (yVar == TRUE) {
        yName <- info$y$name
        xyLoadParTable <- data.frame(
            lhs = paste("t", xName, sep = "_"),
            op = "=~",
            rhs = paste("a", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "d_b"
        )
        yxLoadParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "=~",
            rhs = paste("a", xName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "d_d"
        )
        yyLoadParTable <- data.frame(
            lhs = paste("t", yName, sep = "_"),
            op = "=~",
            rhs = paste("a", yName, 1, sep = "_"),
            user = 1,
            block = 1,
            group = 1,
            free = 1,
            ustart = NA,
            exo = 0,
            label = "d_c"
        )
        loadTable <- do.call(
            rbind,
            list(
                loadParTable,
                xyLoadParTable,
                yxLoadParTable,
                yyLoadParTable
            )
        )
    } else {
        loadTable <- loadParTable
    }
    
    
    loadTable$from <- "dpmLoadings"
    return(loadTable)
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
                        lags = 1,
                        trait = TRUE,
                        state = TRUE,
                        stability = TRUE,
                        crossLag = TRUE,
                        traitCors = TRUE,
                        arCors = TRUE,
                        stateCors = TRUE,
                        residCors = TRUE,
                        predetermined = FALSE,
                        dpm_c = FALSE,
                        dpm_p = FALSE,
                        gclm = FALSE,
                        ma = FALSE,
                        clma = FALSE,
                        slope = "none",
                        measurement = FALSE) {

    ## Check if bivariate or univariate
    yVar <- info$gen$yVar

    if (dpm_c == TRUE | dpm_p == TRUE) {
        trait <- FALSE
    }

    if (gclm == TRUE) {
        trait <- FALSE
    }

    if (trait == FALSE & gclm == FALSE) {
        traitCors == FALSE
    }

    if (ar == FALSE) {
        arCors <- FALSE
    }

    if (state == FALSE) {
        stateCors <- FALSE
    }

    if (dpm_c == TRUE | dpm_p == TRUE | gclm == TRUE) {
        predetermined == FALSE
    }
    
    components <- list(
        ph = .buildPhantom(info),
        obs = .buildObserved(info),
        ar = .buildAr(info),
        imp = .buildImpulses(info),
        residVar = .buildResidVar(info),
        trait = .buildTrait(info),
        dpmTrait = .buildDpmTrait(info),
        gclmTrait = .buildGclmTrait(info),
        stability = .buildStability(info, lags),
        cl = .buildCrossLag(info, lags),
        ma = .buildMa(info),
        clma = .buildClMa(info),
        state = .buildState(info),
        cors = .buildCors(info,
                          ar = arCors,
                          trait = traitCors,
                          state = stateCors,
                          slope = slope),
        preDetCors  = .buildPredetermined(info, slope),
        dpmTraitCors = .buildDpmTraitCors(info),
        dpmCors = .buildDpmCors(info, slope = slope),
        dpmLoadings = .buildDpmLoadings(info),
        residCors = .buildResidCors(info),
        slopes = .buildSlope(info, slope = slope),
        measurement = .buildMeasurement(info)
    )
    
    if (ar == FALSE) {
        components$ar <- NULL
        components$imp <- NULL
    }
    ## Fix resid var stuff here
    if (trait == FALSE) {
        components$trait <- NULL
    }
    if (dpm_c == FALSE & dpm_p == FALSE) {
        components$dpmTrait <- NULL
        components$dpmTraitCors <- NULL
    }
    if (gclm == FALSE) {
        components$gclmTrait <- NULL
    }
    if (stability == FALSE) {
        components$stability <- NULL
    }
    if (ar == FALSE | measurement == TRUE) {
        components$cl <- NULL
    }
    if (ma == FALSE) {
        components$ma <- NULL
    }
    if (clma == FALSE) {
        components$clma <- NULL
    }
    if (state == FALSE) {
        components$state <- NULL
    }
    if (residCors == FALSE) {
        components$residCors <- NULL
    }
    if (predetermined == FALSE) {
        components$preDetCors <- NULL
    }
    if (dpm_p == FALSE) {
        components$dpmCors <- NULL
    }
    if (dpm_c == FALSE) {
        components$dpmLoadings <- NULL
    }
    if (slope == "none") {
        components$slopes <- NULL
    }
    if (measurement == FALSE) {
        components$measurement <- NULL
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
    varNamesInfo$latent <- ifelse(
        varNamesInfo$varName2 == info$x$name &
            (varNamesInfo$varName1 == "l" |
                varNamesInfo$varName1 == "a" |
                varNamesInfo$varName1 == "t" |
                varNamesInfo$varName1 == "s" |
                varNamesInfo$varName1 == "i" |
                varNamesInfo$varName1 == "sl"), 1, NA
    )

    if (yVar == TRUE) {
        varNamesInfo$latent <- ifelse(
        (varNamesInfo$varName2 == info$y$name |
         varNamesInfo$varName2 == info$x$name) &
        (varNamesInfo$varName1 == "l" |
         varNamesInfo$varName1 == "a" |
         varNamesInfo$varName1 == "t" |
         varNamesInfo$varName1 == "s" |
         varNamesInfo$varName1 == "i" |
         varNamesInfo$varName1 == "sl"), 1, NA
        )
    }
        
    latentVarInfo <- varNamesInfo[which(varNamesInfo$latent == 1), 1:4]
    names(latentVarInfo) <- c("varName", "component", "variable", "wave")
    latentVarInfo$indicator <- NA
    if (info$gen$yVar == TRUE) {
        observedVarInfo <- varNamesInfo[which(varNamesInfo$varName1 == info$x$name |
            varNamesInfo$varName1 == info$y$name), 1:4]
    } else {
        observedVarInfo <- varNamesInfo[which(varNamesInfo$varName1 == info$x$name), 1:4]
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
    if (slope == "none") {
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
    } else {
        int.ov <- data.frame(
            lhs = observedVarInfo$varName,
            op = rep("~1", length(observedVarInfo$varName)),
            rhs = rep("", length(observedVarInfo$varName)),
            user = rep(1, length(observedVarInfo$varName)),
            block = rep(1, length(observedVarInfo$varName)),
            group = rep(1, length(observedVarInfo$varName)),
            free = rep(0, length(observedVarInfo$varName)),
            ustart = rep(0, length(observedVarInfo$varName)),
            exo = rep(0, length(observedVarInfo$varName)),
            label = rep("", length(observedVarInfo$varName)),
            from = rep("int.ov", length(observedVarInfo$varName))
        )
    }


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

    if (slope != "none") {
        x_trait_name <- paste("t", info$x$name, sep = "_")
        x_slope_name <- paste("sl", info$x$name, sep = "_")
        int.lv[which(int.lv$lhs == x_trait_name), "free"] <- 1
        int.lv[which(int.lv$lhs == x_slope_name), "free"] <- 1
        int.lv[which(int.lv$lhs == x_trait_name), "ustart"] <- NA
        int.lv[which(int.lv$lhs == x_slope_name), "ustart"] <- NA
        if (info$gen$yVar == TRUE) {
            y_trait_name <- paste("t", info$y$name, sep = "_")
            y_slope_name <- paste("sl", info$y$name, sep = "_")
            int.lv[which(int.lv$lhs == y_trait_name), "free"] <- 1
            int.lv[which(int.lv$lhs == y_slope_name), "free"] <- 1
            int.lv[which(int.lv$lhs == y_trait_name), "ustart"] <- NA
            int.lv[which(int.lv$lhs == y_slope_name), "ustart"] <- NA
        }
    }
    
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
