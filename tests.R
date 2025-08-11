################################################################################
## Create test data
################################################################################

waves <- 10
indicators <- 3
N <- 10000

data_starts_cl <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
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

## Create data frame with correctly labeled variables
data_starts_cl <- data_starts_cl
names(data_starts_cl) <- paste(rep(c("x", "LS"), each = waves),
                      rep(1:waves, 2),
                      sep="_")
for (i in names(data_starts_cl)) {
    data_starts_cl <- addIndicators(data_starts_cl, i, indicators)
}

startC <- (2*waves+1)
endC <- (2*waves)+(2*waves*indicators)
data_starts_cl_I <- data_starts_cl[,startC:endC]
data_starts_cl <- data_starts_cl[,1:(2*waves)]

write_csv(data_starts_cl_I, "tests/data_starts_cl_I.csv")
write_csv(data_starts_cl, "tests/data_starts_cl.csv")

dataI <- read_csv("tests/data_starts_cl_I.csv")
data <- read_csv("tests/data_starts_cl.csv")

m_starts <- panelcoder(data = data,
                       panelModel = "starts",
                       program = "mplus",
                       stationarity = "full")

m_arts <- panelcoder(data = data,
                     panelModel = "arts",
                     program = "mplus",
                     stationarity = "full")


data_starts <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = 0, # Cross lag (Y regressed on X)
    xy = 0, # Cross lag (X regressed on Y)
    cor_xy = .5, # Correlation between X and Y (as correlation)
    xr = 1, # Measurement error for X
    yr = 1 # Measurement error for Y
)

## Create data frame with correctly labeled variables
data_starts <- data_starts
names(data_starts) <- paste(rep(c("x", "LS"), each = waves),
                      rep(1:waves, 2),
                      sep="_")
for (i in names(data_starts)) {
    data_starts <- addIndicators(data_starts, i, indicators)
}

startC <- (2*waves+1)
endC <- (2*waves)+(2*waves*indicators)
data_starts_I <- data_starts[,startC:endC]
data_starts <- data_starts[,1:(2*waves)]

write_csv(data_starts_I, "tests/data_starts_I.csv")
write_csv(data_starts, "tests/data_starts.csv")

dataI <- read_csv("tests/data_starts_I.csv")
data <- read_csv("tests/data_starts.csv")

m_starts <- panelcoder(data = data,
                       panelModel = "starts",
                       program = "mplus",
                       stationarity = "paths")

m_arts <- panelcoder(data = data,
                     panelModel = "arts",
                     program = "mplus",
                     stationarity = "full")
