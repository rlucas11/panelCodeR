getInfo <- function(df) {
    ## Function to get information about names, waves, and indicators from dataframe
    ## Returns errors if data are not structured and named correctly
    ## Warns if length of names exceeds 8 characters (problem for Mplus)
    dfNames <- names(df)

    ## Check length of names for Mplus
    if (max(sapply(tempList, nchar) > 8)) {
        warning("Variable names are longer than 8 characters. This could be a problem if using Mplus")
    }
    
    namesList <- strsplit(dfNames, split = "_")
    nSplits <- sapply(namesList, length)

    ## Check whether splits make sense
    if (any(nSplits > 3)) {
        stop("Too may splits. Check variable names")
    }
    if (length(unique(sapply(namesList, length))) > 1) {
        stop("Variable names cannot be split in the same way for each variable. Check variable names.")
    }
    if (any(nSplits == 1)) {
        stop("Not able to identify number of waves. Check variable names.")
    }
    variableNames <- unique(sapply(namesList, "[[",1))
    if (length(variableNames) == 1) {
        yVar <- FALSE
    } else if (length(variableNames) == 2) {
        yVar <- TRUE
    } else if (length(variableNames) > 2) {
        stop("Variables are not named consistently or there are more than two variables. Check variable names.")
    }

    ## Collect basic info
    xNames <- namesList[sapply(namesList, "[[",1) == variableNames[[1]]]
    xWaves <- as.numeric(unique(sapply(xNames, "[[", 2)))
    if (yVar == TRUE) {
        yNames <- namesList[sapply(namesList, "[[", 1) == variableNames[[2]]]
        yWaves <- as.numeric(unique(sapply(yNames, "[[", 2)))
    }

    ## Check number of indicators
    if (max(nSplits) == 3) {
        if (length(unique(table(sapply(xNames, "[[", 2)))) > 1) {
            stop("Different number of indicators per wave. Check data.")
        }
        xIndicators <- as.numeric(max(unique(sapply(xNames, "[[", 3))))
        if (yVar == TRUE) {
            if (length(unique(table(sapply(yNames, "[[", 2)))) > 1) {
                stop("Different number of indicators per wave. Check data.")
            }
            yIndicators <- as.numeric(max(unique(sapply(yNames, "[[", 3))))
        }
    } else {
        xIndicators <- 1
        if (yVar == TRUE) {
            yIndicators <- 1
        }
    }

    ## Create list of info
    if (yVar == TRUE) {
        info <- list(
            x = list(
                name = variableNames[1],
                waves = xWaves,
                indicators = xIndicators
            ),
            y = list(
                name = variableNames[2],
                waves = yWaves,
                indicators = yIndicators
            )
        )
    } else {
        info <- list(
            x = list(
                name = variableNames[1],
                waves = xWaves,
                indicators = xIndicators
            ),
            y = NULL
        )
    }
    return(info)
}


buildObserved <- function(info) {
    ## Check if bivariate or univariate
    if (is.null(info$y)) {
        yVar <- FALSE
    } else {
        yVar <- TRUE
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

        for (w in info$x$waves) {
            for (i in 1:info$x$indicators) {
                loadingParTable <- data.frame(
                    lhs = paste0("lx", w),
                    op = "=~",
                    rhs = paste("x", w, i, sep="_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = indList$free[i],
                    ustart = indList$ustart[i],
                    exo = 0,
                    label = ""
                )
                var1ParTable <- data.frame(
                    lhs = paste("x", w, i, sep="_"),
                    op = "~~",
                    rhs = paste("x", w, i, sep="_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = ""
                )
                var2ParTable <- data.frame(
                    lhs = paste0("lx", w),
                    op = "=~",
                    rhs = paste0("lx", w),
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
                    var1ParTable,
                    var2ParTable
                )
            }
        }
    } else {
        for (w in info$x$waves) {
            loadingParTable <- data.frame(
                    lhs = paste0("lx", w),
                    op = "=~",
                    rhs = paste("x", w, sep="_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 0,
                    ustart = 1,
                    exo = 0,
                    label = ""
                )
                var1ParTable <- data.frame(
                    lhs = paste("x", w, sep="_"),
                    op = "~~",
                    rhs = paste("x", w, sep="_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = ""
                )
                var2ParTable <- data.frame(
                    lhs = paste0("lx", w),
                    op = "=~",
                    rhs = paste0("lx", w),
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
                    var1ParTable,
                    var2ParTable
                )

            newParTable <- data.frame(
                lhs = paste0("lx", w),
                op = "=~",
                rhs = paste("x", w, sep = "_"),
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
                newParTable
            )
        }
    }

    if (yVar == TRUE) {
        if (info$y$indicators > 1) {
            indList <- list(
                ustart = c(1, rep(NA, (info$y$indicators - 1))),
                free = c(0, rep(1, (info$y$indicators - 1)))
            )

            for (w in info$y$waves) {
                for (i in 1:info$y$indicators) {
                    loadingParTable <- data.frame(
                        lhs = paste0("ly", w),
                        op = "=~",
                        rhs = paste("y", w, i, sep = "_"),
                        user = 1,
                        block = 1,
                        group = 1,
                        free = indList$free[i],
                        ustart = indList$ustart[i],
                        exo = 0,
                        label = ""
                    )
                    var1ParTable <- data.frame(
                        lhs = paste("y", w, i, sep = "_"),
                        op = "~~",
                        rhs = paste("y", w, i, sep = "_"),
                        user = 1,
                        block = 1,
                        group = 1,
                        free = 1,
                        ustart = NA,
                        exo = 0,
                        label = ""
                    )
                    var2ParTable <- data.frame(
                        lhs = paste0("ly", w),
                        op = "=~",
                        rhs = paste0("ly", w),
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
                        var1ParTable,
                        var2ParTable
                    )
                }
            }
        } else {
            for (w in info$y$waves) {
                loadingParTable <- data.frame(
                    lhs = paste0("ly", w),
                    op = "=~",
                    rhs = paste("y", w, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 0,
                    ustart = 1,
                    exo = 0,
                    label = ""
                )
                var1ParTable <- data.frame(
                    lhs = paste("y", w, sep = "_"),
                    op = "~~",
                    rhs = paste("y", w, sep = "_"),
                    user = 1,
                    block = 1,
                    group = 1,
                    free = 1,
                    ustart = NA,
                    exo = 0,
                    label = ""
                )
                var2ParTable <- data.frame(
                    lhs = paste0("ly", w),
                    op = "=~",
                    rhs = paste0("ly", w),
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
                    var1ParTable,
                    var2ParTable
                )

                newParTable <- data.frame(
                    lhs = paste0("ly", w),
                    op = "=~",
                    rhs = paste("y", w, sep = "_"),
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
                    newParTable
                )
            }
        }
    }
    
    return(initialParTable)
}

    ## Create X variable indicators
    
    

