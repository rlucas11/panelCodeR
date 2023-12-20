library(tidyverse)

source("~/Projects/code-generator/scripts/gen_starts.R")
data <- read_csv("testData.csv")

addIndicators <- function(df, var, indicators) {
    var <- sym(var)
    for (i in 1:indicators) {
        label <- letters[i]
        var <- enquo(var)
        prefix <- as_label(var)
        df <- df %>%
            rowwise() %>%
            mutate("{ prefix }{label}" := !!var + rnorm(1, 0, 1))
    }
    return(df)
}


for (i in names(data)) {
    data <- addIndicators(data, i, 3)
}

test <- run_starts_mplus(data[1:1000,],
                         5,
#                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 1
                         )

test <- run_arts_mplus(data[1:1000,],
                         5,
                         xWaves = c(1:5, 7:10),
                         xIndicators = 3,
                         yIndicators = 1
                         )


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

