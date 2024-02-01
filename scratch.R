################################################################################
## Setup
################################################################################

library(tidyverse)
library(lavaan)
library(MplusAutomation)
load_all()
#library(panelCodeR)

data <- gen_starts(
    n = 1000, # N to generate
    nwaves = 10, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = .2, # Cross lag (Y regressed on X)
    xy = .1, # Cross lag (X regressed on Y)
    cor_xy = .5, # Correlation between X and Y (as correlation)
    xr = 1, # Measurement error for X
    yr = 1 # Measurement error for Y
)



data1 <- gen_starts(
    n = 1000, # N to generate
    nwaves = 2, # Number of waves
    ri_x = 0, # Random intercept variance for X
    ri_y = 0, # Random intercept variance for Y
    cor_i = .5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = -.2, # Cross lag (Y regressed on X)
    xy = -.2, # Cross lag (X regressed on Y)
    cor_xy = -.5, # Correlation between X and Y (as correlation)
    xr = .2, # Measurement error for X
    yr = .2 # Measurement error for Y
)


summary(lm(y2 ~ y1 + x1, data=data1))
summary(lm(y1 ~ x1 + y2, data=data1))


data2 <- gen_starts(
    n = 10000, # N to generate
    nwaves = 2, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = -.5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = 0, # Cross lag (Y regressed on X)
    xy = 0, # Cross lag (X regressed on Y)
    cor_xy = -.5, # Correlation between X and Y (as correlation)
    xr = .2, # Measurement error for X
    yr = .2 # Measurement error for Y
)

summary(lm(y2 ~ y1 + x1, data=data2))
summary(lm(y1 ~ x1 + y2, data=data2))


data3 <- gen_starts(
    n = 10000, # N to generate
    nwaves = 2, # Number of waves
    ri_x = 0, # Random intercept variance for X
    ri_y = 0, # Random intercept variance for Y
    cor_i = -.5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = 0, # Cross lag (Y regressed on X)
    xy = 0, # Cross lag (X regressed on Y)
    cor_xy = -.5, # Correlation between X and Y (as correlation)
    xr = .2, # Measurement error for X
    yr = .2 # Measurement error for Y
)

summary(lm(y2 ~ y1 + x1, data=data2))
summary(lm(y1 ~ x1 + y2, data=data2))



################################################################################
## Testing Rewrite
################################################################################


startsModel <- buildLavaan(waves = 5)
oldParTable <- lavaanify(startsModel)



## Create data frame with correctly labeled variables
data2 <- data
names(data2) <- paste(rep(c("X", "Y"), each = 10),
                      rep(1:10, 2),
                      sep="_")
for (i in names(data2)) {
    data2 <- addIndicators(data2, i, 3)
}
dataI <- data2[,21:80]
data2 <- data2[,1:20]



## Testing
infoI <- getInfo(dataI)
info <- getInfo(data2)


out <- panelcoder(data=data2)

##out <- panelcoder(data=data2[,c(1:5,11:15)], panelModel="starts")
summary(out[[2]])

lavExport(out, target="mplus")

outM <- run_starts_mplus(data, waves=10)

testM <- lav2mplus(out[[1]])

mplusStatement <- mplusObject(TITLE="test",
                              rdata=data2,
                              ANALYSIS="MODEL=NOCOVARIANCES;",
                              MODEL=testM)

mplusModeler(mplusStatement, modelout = "testPc.inp", run=1)


## Phantom vars

panelcoder(data2[,c(1:5,11,12,14,15)])


## Latent vars
info <- getInfo(dataI)

latent <- panelcoder(dataI[,c(1:15,31:45)])
testL <- lav2mplus(latent[[1]])

mplusStatement <- mplusObject(TITLE="test",
                              rdata=dataI[,c(1:15, 31:45)],
                              ANALYSIS="MODEL=NOCOVARIANCES;",
                              MODEL=testL)

mplusModeler(mplusStatement, modelout = "testPcL.inp", run=1)


## Mix of latent and observed

testMix <- panelcoder(data2[,c(1:5, 51:65)])
testMixM <- lav2mplus(testMix[[1]])

mplusStatement <- mplusObject(TITLE="test",
                              rdata=data2[,c(1:5, 51:65)],
                              ANALYSIS="MODEL=NOCOVARIANCES;",
                              MODEL=testMixM)

mplusModeler(mplusStatement, modelout = "testPcL.inp", run=1)


################################################################################
################################################################################


testL <- run_starts_lavaan(data,5)

testL <- run_startsY_lavaan(data[1:1000,],5)

test <- run_dpm_mplus(data,
                      10
                      )



test <- buildMplus(
                   waves = 5
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


testM <- run_starts_mplus(data,
                         5
                         )
summary(testM)

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

testLavaan <- run_starts_lavaan(data = data, waves = 10)

## Bivariate STARTS

startsModel <- buildLavaan(waves = 5,
                           xWaves = 1:5,
                           yWaves = 1:5,
                           xIndicators = 1,
                           yIndicators = 1,
                           stateCor = FALSE
                           )
cat(startsModel)


startsModel <- buildLavaan(waves = 5)
cat(startsModel)
startsFit <- lavaan(startsModel, data)
summary(startsFit)                             

startsFit2 <- lavaan(startsModel, data, missing="FIML", meanstructure=TRUE, int.lv.free=FALSE)
summary(startsFit2)

run_starts_mplus(data, waves=5)

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

startsModelY <- lavaanStartsY(5, 1:5)
cat(startsModelY)
startsFitY <- lavaan(startsModelY, data[1:1000,])
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
                      10
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


testDpm <- run_dpm_lavaan(data,
                          waves = 5)


################################################################################
## rewrite
################################################################################

startsModel <- buildLavaan(waves = 5)
oldParTable <- lavaanify(startsModel)

## Create data frame with correctly labeled variables
data2 <- data
names(data2) <- paste(rep(c("X", "Y"), each = 10),
                      rep(1:10, 2),
                      sep="_")
info <- getInfo(data2[1:1000, c(1:6, 11:16)])


for (i in names(data2)) {
    data2 <- addIndicators(data2, i, 3)
}

dataI <- data2[,21:80]
data2 <- data2[,1:20]



## Testing
infoI <- getInfo(dataI)
info <- getInfo(data2)


panelcoder(data=data2)


model <- .buildTable(info, stateCors = FALSE)
lavModel <- model$model

lavModel <- .constrainStability(model$model, info)
lavModel <- .constrainCl(lavModel, info)
lavModel <- .arStationarity(lavModel, info)
lavModel <- .constrainStateVar(lavModel, info)
lavModel <- .constrainStateCors(lavModel, info, zero = TRUE)
lavModel <- .buildLimits(lavModel, info)


temp <- lavaan(model = lavModel, data = data2)
summary(temp)

model <- .buildTable(infoI, stateCors = FALSE)
lavModel <- model$model

lavModel <- .constrainStability(lavModel, info)
lavModel <- .constrainCl(lavModel, info)
lavModel <- .arStationarity(lavModel, info)
lavModel <- .constrainStateVar(lavModel, info)
lavModel <- .constrainStateCors(lavModel, info, zero = TRUE)
lavModel <- .buildLimits(lavModel, info, stateCors = FALSE)


temp <- lavaan::lavaan(model = lavModel, data = dataI)
summary(temp)


################################################################################
## Test simple model
################################################################################

simpleData <- data2[1:500, c(1:5, 11:15)]
simpleDataO <- data[1:500, c(1:5, 11:15)]

simpleInfo <- getInfo(simpleData)


simple <- .buildTable(simpleInfo, trait = TRUE, state = TRUE, stateCors = FALSE)
simple <- lavaan(simple$model, simpleData)

out <- run_riclpm_lavaan(simpleDataO, 3, stationarity = FALSE)

lavExport(simple, target = "mplus", prefix = "test")


simple2 <- panelcoder(simpleData, "starts")
lavExport(simple2, target="mplus", prefix="test2")


testTable <- .buildConstraint(simple$model,
                              "a2",
                              "==",
                              "a3")
testTable <- .buildConstraint(testTable,
                              "b2",
                              "==",
                              "b3")
testTable <- .buildConstraint(testTable,
                              "c2",
                              "==",
                              "c3")
testTable <- .buildConstraint(testTable,
                              "d2",
                              "==",
                              "d3")


testTable2 <- .buildConstraint(testTable,
                               "xvar2",
                               "==",
                               "xvar1 - a2*a2*xvar1 - d2*d2*yvar1 - 2*a2*d2*cov_ar1")
testTable2 <- .buildConstraint(testTable2,
                               "yvar2",
                               "==",
                               "yvar1 - b2*b2*yvar1 - c2*c2*xvar1 - 2*b2*c2*cov_ar1")
testTable2 <- .buildConstraint(testTable2,
                              "xvar2",
                              "==",
                              "xvar3")
testTable2 <- .buildConstraint(testTable2,
                              "yvar2",
                              "==",
                              "yvar3")





simple <- .buildTable(simpleInfo, trait = TRUE, state = TRUE, stateCors = FALSE)

testTable <- .constrainStability(simple$model, simpleInfo)
testTable <- .constrainCl(testTable, simpleInfo)
testTable <- .arStationarity(testTable, simpleInfo)
testTable <- .buildLimits(testTable, simpleInfo)
testTable <- .constrainStateVar(testTable, simpleInfo)

test <- lavaan(testTable, simpleData)
summary(test)

lavExport(test, target="mplus")

test2 <- run_riclpm_lavaan(simpleDataO, 5)
test2 <- run_starts_lavaan(simpleDataO, 5)

test3 <- run_starts_mplus(simpleDataO, 5)



################################################################################
## Scratch
################################################################################


testStarts <- buildLavaan(waves=3, YVar=FALSE, trait=FALSE, state=FALSE)

test <- lavaan(testModel, data)
summary(test)

test <- lavaan(testStarts, data, meanstructure=TRUE, int.ov.free=TRUE, int.lv.free=FALSE)
summary(test)



testModel <- "
x2 ~ x1
"

test <- lavaan(testModel,
               data,
               meanstructure=TRUE,
               int.ov.free=TRUE,    ## Intercepts of observed free
               int.lv.free=FALSE,    ## Intercepts of latent free
               auto.fix.first=TRUE,
               auto.fix.single=TRUE,
               auto.var=TRUE,
               auto.cov.lv.x=TRUE,
               auto.efa=TRUE,
               auto.th=TRUE,
               auto.delta=TRUE,
               auto.cov.y=TRUE)

summary(test)


components <- list(
    obs = .buildObserved(info),
    ar = .buildAr(info),
    trait = .buildTrait(info),
    stability = .buildStability(info),
    cl = .buildCrossLag(info),
    state = .buildState(info),
    cors = .buildCors(info, ar = ar, trait = trait, state = stateCors),
    residCors = .buildResidCors(info)
)
