getInfo <- function(df) {
    ## Function to get information about names, waves, and indicators from dataframe
    ## Returns errors if data are not structured and named correctly
    ## Warns if length of names exceeds 8 characters (problem for Mplus)
    dfNames <- names(df)

    ## Check length of names for Mplus
    if (max(sapply(dfNames, nchar) > 8)) {
        warning("Variable names are longer than 8 characters. This could be a problem if using Mplus",
                call. = FALSE)
    }

    ## Collect basic info
    namesList <- strsplit(dfNames, split = "_")
    nSplits <- sapply(namesList, length)

    variableNames <- unique(sapply(namesList, "[[",1))
    if (length(variableNames) == 1) {
        yVar <- FALSE
    } else if (length(variableNames) == 2) {
        yVar <- TRUE
    } else if (length(variableNames) > 2) {
        stop("Variables are not named consistently or there are more than two variables. Check variable names.",
             call. = FALSE)
    }


    xNames <- namesList[sapply(namesList, "[[",1) == variableNames[[1]]]
    xWaves <- as.numeric(unique(sapply(xNames, "[[", 2)))
    xSplits <- sapply(xNames, length)
    maxWaves <- max(xWaves)
    if (yVar == TRUE) {
        yNames <- namesList[sapply(namesList, "[[", 1) == variableNames[[2]]]
        ySplits <- sapply(yNames, length)
        yWaves <- as.numeric(unique(sapply(yNames, "[[", 2)))
        maxWaves <- max(xWaves, yWaves)
    }


    ## Check whether splits make sense
    if (any(nSplits > 3)) {
        stop("Too may splits. Check variable names", call. = FALSE)
    }
    if (length(unique(sapply(xNames, length))) > 1) {
        stop("X variable names cannot be split in the same way for each variable. Check variable names.",
             call. = FALSE)
    }
    if (yVar == TRUE) {
        if (length(unique(sapply(yNames, length))) > 1) {
            stop("Y variable names cannot be split in the same way for each variable. Check variable names.",
                 call. = FALSE)
        }
    }
    if (any(nSplits == 1)) {
        stop("Not able to identify number of waves. Check variable names.", call.=FALSE)
    }
    
    

    ## Check number of indicators
    if (max(xSplits) == 3) {
        if (length(unique(table(sapply(xNames, "[[", 2)))) > 1) {
            stop("Different number of indicators per wave. Check data.", call. = FALSE)
        }
        xIndicators <- as.numeric(max(unique(sapply(xNames, "[[", 3))))
    } else {
        xIndicators <- 1
    }
    
    if (yVar == TRUE) {
        if (max(ySplits) == 3) {
            if (length(unique(table(sapply(yNames, "[[", 2)))) > 1) {
                stop("Different number of indicators per wave. Check data.", call. = FALSE)
            }
            yIndicators <- as.numeric(max(unique(sapply(yNames, "[[", 3))))
        } else {
            yIndicators <- 1
        }
    } 

    ## Create list of info
    if (yVar == TRUE) {
        info <- list(
            gen = list(
                maxWaves = maxWaves,
                yVar = TRUE
                ),
            x = list(
                name = variableNames[1],
                waves = 1:maxWaves,
                actualWaves = xWaves,
                indicators = xIndicators
            ),
            y = list(
                name = variableNames[2],
                waves = 1:maxWaves,
                actualWaves = yWaves,
                indicators = yIndicators
            )
        )
    } else {
        info <- list(
            gen = list(
                maxWaves = maxWaves,
                yVar = FALSE
                ),
            x = list(
                name = variableNames[1],
                waves = xWaves,
                actualWaves = xWaves,
                indicators = xIndicators
            ),
            y = NULL
        )
    }
    return(info)
}


