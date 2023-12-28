################################################################################
## Setup
################################################################################

library(tidyverse)
library(lavaan)
library(panelCodeR)

data <- gen_starts(
    n = 50000, # N to generate
    nwaves = 10, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = .1, # Cross lag (Y regressed on X)
    xy = .05, # Cross lag (X regressed on Y)
    cor_xy = .5, # Correlation between X and Y (as correlation)
    xr = 1, # Measurement error for X
    yr = 1 # Measurement error for Y
)


test <- buildMplus(data,
                   waves = 5,
                   xWaves = c(1,2,4,5)
                   )

test <- run_starts_mplus(data,
                         5,
                         YVar = FALSE,
                         trait = FALSE,
                         state = FALSE
                         )

summary(test)



## Load sample data
data <- read_csv("testData.csv")

for (i in names(data)) {
    data <- addIndicators(data, i, 3)
}

################################################################################
## Mplus
################################################################################


test <- run_starts_mplus(data,
                         5,
                         stateCor = TRUE
                         )
summary(test)

test <- run_starts_mplus(data,
                         5,
                         xIndicators = 3,
                         yIndicators = 3,
                         stateCor = TRUE)
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
test$results$parameters$stdyx.standardized

test <- run_riclpm_mplus(data,
                         5
                         )

summary(test)
test$results$parameters$unstandardized


test <- run_clpm_mplus(data,
                       5
                       )


summary(test)
test$results$parameters$unstandardized

test <- run_startsx_mplus(data,
                         10
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

## Bivariate STARTS

startsModel <- buildLavaan(waves = 5,
                           xWaves = 1:5,
                           yWaves = 1:5,
                           xIndicators = 3,
                           yIndicators = 3,
                           stateCor = TRUE
                           )
cat(startsModel)

startsFit <- lavaan(startsModel, data)
summary(startsFit)                             

startsModel <- lavaanStarts2(5, 1:5, 1:5, stateCor = FALSE)
startsFit <- lavaan(startsModel, data)
summary(startsFit)

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
startsFit2 <- lavaan(startsMod2, data, missing = "FIML",
    estimator = "MLR")
summary(startsFit2)


riclpmMod <- lavaanRiclpm(5, 1:5, 1:5)
cat(riclpmMod)
riclpmFit <- lavaan(riclpmMod, data)
summary(riclpmFit)

clpmMod <- lavaanClpm(5, 1:5, 1:5)
cat(clpmMod)
clpmFit <- lavaan(clpmMod, data)
summary(clpmFit)


artsMod <- lavaanArts(10, 1:10)
cat(artsMod)
artsFit <- lavaan(artsMod, data)
summary(artsFit)


################################################################################
## DPM
################################################################################

test <- run_dpm_mplus(data,
                      5
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

################################################################################
## work on class
################################################################################


outputList <- list(mplusOutput = output,
    xTrait = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "TRAIT.X"),
        "est"
    ],
    xAr = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "AR.P.X"),
        "est"
    ],
    xState = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "ST.P.X"),
        "est"
    ],
    xStability = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "STAB.X"),
        "est"
    ],
    yTrait = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "TRAIT.Y"),
        "est"
    ],
    yAr = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "AR.P.Y"),
        "est"
    ],
    yState = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "ST.P.Y"),
        "est"
    ],
    yStability = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "STAB.Y"),
        "est"
    ],
    traitCor = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "COR_TXTY"),
        "est"
    ],
    arCor = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "COR_ARXA"),
        "est"
    ],
    stateCor = output$results$parameters$unstandardized[
        which(output$results$parameters$unstandardized$param == "COR_SXSY"),
        "est"
    ]
)

class(testListS3) <- "pcObject"

summary.pcObject <- function(obj) {
    cat("Model Summary: \n")
    summary(obj$mplusOutput)
    cat("\n")
    cat("Variance Decomposition for X Variable: \n")
    cat("\n")
    cat("Trait: ", obj$xTrait,
        ", Autoregressive: ",obj$xAr,
        ", State: ",obj$xState, "\n")
    cat("\n")
    cat("Stability of X: ", obj$xStab, "\n")
    cat("\n")
    cat("Variance Decomposition for Y Variable: \n")
    cat("\n")
    cat("Trait: ", obj$yTrait,
        ", Autoregressive: ", obj$yAr,
        ", State: ", obj$yState, "\n")
    cat("\n")
    cat("Stability of Y: ", obj$yStab, "\n")
    cat("\n")
    cat("Correlations: \n")
    cat("\n")
    cat("Stable Trait: ", obj$traitCor, "\n")
    cat("Autoregressive Trait: ", obj$arCor, "\n")
    cat("State: ", obj$stateCor, "\n")
    }
    
################################################################################
## compareUnivariate
################################################################################

compareUnivariate(data, 10)

################################################################################
## Lavaan DPM
################################################################################

dpmModel <- buildLavaanDpm(waves = 5)
cat(dpmModel)

dpmFit <- lavaan(dpmModel, data)
summary(dpmFit)                             
