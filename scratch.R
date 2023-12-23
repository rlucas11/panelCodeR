################################################################################
## Setup
################################################################################

library(tidyverse)
library(lavaan)
library(panelCodeR)

## For creating new data
source("~/Projects/code-generator/scripts/gen_starts.R")

tempData <- gen_starts(
    n = 500, # N to generate
    nwaves = 10, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = .4, # Cross lag (Y regressed on X)
    xy = .2, # Cross lag (X regressed on Y)
    cor_xy = .5, # Correlation between X and Y (as correlation)
    xr = 0, # Measurement error for X
    yr = 0 # Measurement error for Y
)



## Load sample data
data <- read_csv("testData.csv")

addIndicators <- function(df, var, indicators) {
    var <- rlang::sym(var)
    for (i in 1:indicators) {
        label <- letters[i]
        var <- rlang::enquo(var)
        prefix <- rlang::as_label(var)
        df <- df %>%
            dplyr::rowwise() %>%
            dplyr::mutate("{ prefix }{label}" := !!var + rnorm(1, 0, 1))
    }
    return(df)
}


for (i in names(data)) {
    data <- addIndicators(data, i, 3)
}

################################################################################
## Mplus
################################################################################


test <- run_starts_mplus(data[1:1000,],
                         5)
summary(test)

test <- run_starts_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 3)
summary(test)



test <- run_starts_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 1,
                         constrainCors = FALSE
                         )
summary(test)

test <- run_arts_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 1
                         )
summary(test)

test <- run_sts_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 2
                         )
summary(test)
test$results$parameters$unstandardized

test <- run_riclpm_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 2
                         )

summary(test)
test$results$parameters$unstandardized


test <- run_clpm_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 2
                         )

summary(test)
test$results$parameters$unstandardized

test <- run_startsx_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3
                         )

test$results$parameters$unstandardized
summary(test)


test <- run_startsx_mplus(data[1:1000,],
                          5,
                          xWaves = c(1,2,3,5),
                          xIndicators = 3
                          )

run_startsy_mplus(data[1:1000,],
                  5,
                  yWaves = c(1,2,3,5),
                  yIndicators = 1
                  )

test$results$parameters$unstandardized
summary(test)

################################################################################
## lavaan
################################################################################

## Univariate Starts
startsModelX <- lavaanStartsX(10, 1:10)
cat(startsModelX)
startsFitX <- lavaan(startsModelX, data)
summary(startsFitX)

startsModelX <- lavaanStartsX(10, xWaves = c(1,2,3,5,6,7,9,10))
cat(startsModelX)
startsFitX <- lavaan(startsModelX, data)
summary(startsFitX)

startsModelY <- lavaanStartsY(10, 1:10)
cat(startsModelY)
startsFitY <- lavaan(startsModelY, data)
summary(startsFitY)
fitMeasures(startsFitY)

startsMod2 <- lavaanStarts2(10, xWaves = c(1,2,3,5,6,7,9,10), yWaves = c(2,3,5,6,8,9))
cat(startsMod2)
startsFit2 <- lavaan(startsMod2, data)
summary(startsFit2)


riclpmMod <- lavaanRiclpm(10, 1:10)
cat(riclpmMod)
riclpmFit <- lavaan(riclpmMod, data[1:1000,])
summary(riclpmFit)

clpmMod <- lavaanClpm(10,1:10)
cat(clpmMod)
clpmFit <- lavaan(clpmMod, data[1:1000,])
summary(clpmFit)


artsMod <- lavaanArts(10, 1:10)
cat(artsMod)
artsFit <- lavaan(artsMod, data)
summary(artsFit)


################################################################################
## DPM
################################################################################

test <- run_dpm_mplus(data[1:1000,],
                      4,
                      xIndicators = 3,
                      constrainCors = FALSE
                      )

summary(test)
test$results$parameters$unstandardized
test$results$parameters$stdyx.standardized


test <- run_dpm_mplus(data[1:1000,],
                      5,
                      xIndicators = 3,
                      yIndicators = 1,
                      state = TRUE,
                      analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001; ITERATIONS=20000"
                      )

summary(test)
test$results$parameters$unstandardized
test$results$parameters$stdyx.standardized
