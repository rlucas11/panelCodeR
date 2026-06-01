################################################################################
## Setup
################################################################################

library(tidyverse)
## library(lavaan)
## library(MplusAutomation)

load_all()

## library(panelCodeR)

dataI <- read_csv("tests/testDataI.csv")
data2 <- read_csv("tests/testData2.csv")
data <- read_csv("tests/test_data.csv")


test1 <- panelcoder(data[1:10],
                   panelModel="starts",
                   program = "mplus",
                   stationarity="paths",
                   constrainState = TRUE,
                   limits = TRUE)

test <- panelcoder(data2[1:500,c(1:7,11:17)],
                   panelModel="arts",
                   program = "mplus",
                   stationarity="paths",
                   constrainState = TRUE)



test <- panelcoder(dataI, panelModel="starts", program = "mplus", stationarity="full")
test <- panelcoder(data2, panelModel="starts", program = "mplus", stationarity = "full")

test <- panelcoder(data2, panelModel="riclpm", program = "mplus", stationarity = "full", title="riclpm")
test <- panelcoder(data2, panelModel="arts", program = "mplus", stationarity = "full", title="arts")
test <- panelcoder(data2, panelModel="dpm_p", program = "mplus", stationarity = "full")
test <- panelcoder(data2, panelModel="dpm_c", program = "mplus", stationarity = "full")

testLavaan <- panelcoder(data2[c(1:7,11:17)],
                         panelModel="starts",
                         program = "lavaan",
                         stationarity = "full")



model <- .buildModel(data2)
test <- lavaan(model, data2, missing='fiml')

test <- panelcoder(data2, panelModel="starts", program="mplus")

test2 <- panelcoder(dataI, panelModel="starts", program="mplus")

test4 <- panelcoder(dataI,
                    panelModel="starts",
                    program = "lavaan",
                    invariance = TRUE,
                    residVar = TRUE)

test2m <- panelcoder(dataI, panelModel="starts", program="lavaan", invariance=TRUE, residCors=TRUE)


testInfo <- getInfo(dataI)


create_custom_matrix <- function(rows, cols) {
  # Initialize an empty matrix
  mat <- matrix(0, nrow = rows, ncol = cols)
  
  # Fill the matrix with sequential values
  counter <- 1
  for (i in 1:rows) {
    for (j in 1:cols) {
      mat[i, j] <- counter
      counter <- counter + 1
    }
  }
  
  # Adjust values in the matrix according to the specified pattern
  for (j in 2:cols) {
    mat[2:rows, j] <- mat[2:rows, j] + (cols - j + 1) * rows
  }
  
  return(mat)
}

create_custom_matrix(4,5)

create_matrix <- function(rows, cols) {
  # Initialize an empty matrix
  mat <- matrix(0, nrow = rows, ncol = cols)
  
  # Fill the matrix with sequential values
  counter <- 1
  for (i in 1:rows) {
    for (j in 1:cols) {
      mat[i, j] <- counter
      counter <- counter + 1
    }
  }
  
  # Reverse the even rows in the matrix
  mat[seq(2, rows, by = 2), ] <- mat[seq(2, rows, by = 2), ncol(mat):1]
  
  return(mat)
}


split_integer <- function(number, parts) {
  quotient <- floor(number / parts)
  remainder <- number %% parts
  
  result <- rep(quotient, times = parts)
  
  if (remainder > 0) {
    result[parts] <- result[parts] + remainder
  }
  
  return(result)
}

# Example usage:
integer_to_split <- 14
num_parts <- 3
result_vector <- split_integer(integer_to_split, num_parts)

# Display the result
print(result_vector)


split_integer <- function(number, parts) {
  quotient <- floor(number / parts)
  remainder <- number %% parts
  
  result <- rep(quotient, times = parts)
  
  if (remainder > 0) {
    result[1:remainder] <- result[1:remainder] + 1
  }
  
  return(result)
}


testModel <- '
l_x_1 =~ 1*x_1
x_1 ~~ 0*x_1
l_x_2 =~ 1*x_2
x_2 ~~ 0*x_2
l_x_3 =~ 1*x_3
x_3 ~~ 0*x_3
l_x_4 =~ 1*x_4
x_4 ~~ 0*x_4
l_x_5 =~ 1*x_5
x_5 ~~ 0*x_5
l_x_6 =~ 1*x_6
x_6 ~~ 0*x_6
l_x_7 =~ 1*x_7
x_7 ~~ 0*x_7
l_x_1 ~~ 0*l_x_1
l_x_2 ~~ 0*l_x_2
l_x_3 ~~ 0*l_x_3
l_x_4 ~~ 0*l_x_4
l_x_5 ~~ 0*l_x_5
l_x_6 ~~ 0*l_x_6
l_x_7 ~~ 0*l_x_7
l_LS_1 =~ 1*LS_1
LS_1 ~~ 0*LS_1
l_LS_2 =~ 1*LS_2
LS_2 ~~ 0*LS_2
l_LS_3 =~ 1*LS_3
LS_3 ~~ 0*LS_3
l_LS_4 =~ 1*LS_4
LS_4 ~~ 0*LS_4
l_LS_5 =~ 1*LS_5
LS_5 ~~ 0*LS_5
l_LS_6 =~ 1*LS_6
LS_6 ~~ 0*LS_6
l_LS_7 =~ 1*LS_7
LS_7 ~~ 0*LS_7
l_LS_1 ~~ 0*l_LS_1
l_LS_2 ~~ 0*l_LS_2
l_LS_3 ~~ 0*l_LS_3
l_LS_4 ~~ 0*l_LS_4
l_LS_5 ~~ 0*l_LS_5
l_LS_6 ~~ 0*l_LS_6
l_LS_7 ~~ 0*l_LS_7
a_x_1 =~ 1*l_x_1
a_x_2 =~ 1*l_x_2
a_x_3 =~ 1*l_x_3
a_x_4 =~ 1*l_x_4
a_x_5 =~ 1*l_x_5
a_x_6 =~ 1*l_x_6
a_x_7 =~ 1*l_x_7
a_LS_1 =~ 1*l_LS_1
a_LS_2 =~ 1*l_LS_2
a_LS_3 =~ 1*l_LS_3
a_LS_4 =~ 1*l_LS_4
a_LS_5 =~ 1*l_LS_5
a_LS_6 =~ 1*l_LS_6
a_LS_7 =~ 1*l_LS_7
a_x_1 ~~ 0*a_x_1
a_x_2 ~~ 0*a_x_2
a_x_3 ~~ 0*a_x_3
a_x_4 ~~ 0*a_x_4
a_x_5 ~~ 0*a_x_5
a_x_6 ~~ 0*a_x_6
a_x_7 ~~ 0*a_x_7
a_LS_1 ~~ 0*a_LS_1
a_LS_2 ~~ 0*a_LS_2
a_LS_3 ~~ 0*a_LS_3
a_LS_4 ~~ 0*a_LS_4
a_LS_5 ~~ 0*a_LS_5
a_LS_6 ~~ 0*a_LS_6
a_LS_7 ~~ 0*a_LS_7
i_x_1 =~ 1*a_x_1
i_x_2 =~ 1*a_x_2
i_x_3 =~ 1*a_x_3
i_x_4 =~ 1*a_x_4
i_x_5 =~ 1*a_x_5
i_x_6 =~ 1*a_x_6
i_x_7 =~ 1*a_x_7
i_LS_1 =~ 1*a_LS_1
i_LS_2 =~ 1*a_LS_2
i_LS_3 =~ 1*a_LS_3
i_LS_4 =~ 1*a_LS_4
i_LS_5 =~ 1*a_LS_5
i_LS_6 =~ 1*a_LS_6
i_LS_7 =~ 1*a_LS_7
i_x_1 ~~ xvar1*i_x_1
i_x_2 ~~ xvar2*i_x_2
i_x_3 ~~ xvar3*i_x_3
i_x_4 ~~ xvar4*i_x_4
i_x_5 ~~ xvar5*i_x_5
i_x_6 ~~ xvar6*i_x_6
i_x_7 ~~ xvar7*i_x_7
i_LS_1 ~~ yvar1*i_LS_1
i_LS_2 ~~ yvar2*i_LS_2
i_LS_3 ~~ yvar3*i_LS_3
i_LS_4 ~~ yvar4*i_LS_4
i_LS_5 ~~ yvar5*i_LS_5
i_LS_6 ~~ yvar6*i_LS_6
i_LS_7 ~~ yvar7*i_LS_7
t_x =~ 1*l_x_1
t_x =~ 1*l_x_2
t_x =~ 1*l_x_3
t_x =~ 1*l_x_4
t_x =~ 1*l_x_5
t_x =~ 1*l_x_6
t_x =~ 1*l_x_7
t_LS =~ 1*l_LS_1
t_LS =~ 1*l_LS_2
t_LS =~ 1*l_LS_3
t_LS =~ 1*l_LS_4
t_LS =~ 1*l_LS_5
t_LS =~ 1*l_LS_6
t_LS =~ 1*l_LS_7
t_x ~~ x_tVar*t_x
t_LS ~~ y_tVar*t_LS
a_x_2 ~ a2*a_x_1
a_x_3 ~ a3*a_x_2
a_x_4 ~ a4*a_x_3
a_x_5 ~ a5*a_x_4
a_x_6 ~ a6*a_x_5
a_x_7 ~ a7*a_x_6
a_LS_2 ~ b2*a_LS_1
a_LS_3 ~ b3*a_LS_2
a_LS_4 ~ b4*a_LS_3
a_LS_5 ~ b5*a_LS_4
a_LS_6 ~ b6*a_LS_5
a_LS_7 ~ b7*a_LS_6
a_x_2 ~ d2*a_LS_1
a_LS_2 ~ c2*a_x_1
a_x_3 ~ d3*a_LS_2
a_LS_3 ~ c3*a_x_2
a_x_4 ~ d4*a_LS_3
a_LS_4 ~ c4*a_x_3
a_x_5 ~ d5*a_LS_4
a_LS_5 ~ c5*a_x_4
a_x_6 ~ d6*a_LS_5
a_LS_6 ~ c6*a_x_5
a_x_7 ~ d7*a_LS_6
a_LS_7 ~ c7*a_x_6
s_x_1 =~ 1*l_x_1
s_x_2 =~ 1*l_x_2
s_x_3 =~ 1*l_x_3
s_x_4 =~ 1*l_x_4
s_x_5 =~ 1*l_x_5
s_x_6 =~ 1*l_x_6
s_x_7 =~ 1*l_x_7
s_LS_1 =~ 1*l_LS_1
s_LS_2 =~ 1*l_LS_2
s_LS_3 =~ 1*l_LS_3
s_LS_4 =~ 1*l_LS_4
s_LS_5 =~ 1*l_LS_5
s_LS_6 =~ 1*l_LS_6
s_LS_7 =~ 1*l_LS_7
s_x_1 ~~ sx1*s_x_1
s_x_2 ~~ sx2*s_x_2
s_x_3 ~~ sx3*s_x_3
s_x_4 ~~ sx4*s_x_4
s_x_5 ~~ sx5*s_x_5
s_x_6 ~~ sx6*s_x_6
s_x_7 ~~ sx7*s_x_7
s_LS_1 ~~ sy1*s_LS_1
s_LS_2 ~~ sy2*s_LS_2
s_LS_3 ~~ sy3*s_LS_3
s_LS_4 ~~ sy4*s_LS_4
s_LS_5 ~~ sy5*s_LS_5
s_LS_6 ~~ sy6*s_LS_6
s_LS_7 ~~ sy7*s_LS_7
i_x_1 ~~ cov_ar1*i_LS_1
i_x_2 ~~ cov_ar2*i_LS_2
i_x_3 ~~ cov_ar3*i_LS_3
i_x_4 ~~ cov_ar4*i_LS_4
i_x_5 ~~ cov_ar5*i_LS_5
i_x_6 ~~ cov_ar6*i_LS_6
i_x_7 ~~ cov_ar7*i_LS_7
t_x ~~ cov_txty*t_LS
x_1 ~ 1
x_2 ~ 1
x_3 ~ 1
x_4 ~ 1
x_5 ~ 1
x_6 ~ 1
x_7 ~ 1
LS_1 ~ 1
LS_2 ~ 1
LS_3 ~ 1
LS_4 ~ 1
LS_5 ~ 1
LS_6 ~ 1
LS_7 ~ 1
l_x_1 ~ 0*1
l_x_2 ~ 0*1
l_x_3 ~ 0*1
l_x_4 ~ 0*1
l_x_5 ~ 0*1
l_x_6 ~ 0*1
l_x_7 ~ 0*1
l_LS_1 ~ 0*1
l_LS_2 ~ 0*1
l_LS_3 ~ 0*1
l_LS_4 ~ 0*1
l_LS_5 ~ 0*1
l_LS_6 ~ 0*1
l_LS_7 ~ 0*1
a_x_1 ~ 0*1
a_x_2 ~ 0*1
a_x_3 ~ 0*1
a_x_4 ~ 0*1
a_x_5 ~ 0*1
a_x_6 ~ 0*1
a_x_7 ~ 0*1
a_LS_1 ~ 0*1
a_LS_2 ~ 0*1
a_LS_3 ~ 0*1
a_LS_4 ~ 0*1
a_LS_5 ~ 0*1
a_LS_6 ~ 0*1
a_LS_7 ~ 0*1
i_x_1 ~ 0*1
i_x_2 ~ 0*1
i_x_3 ~ 0*1
i_x_4 ~ 0*1
i_x_5 ~ 0*1
i_x_6 ~ 0*1
i_x_7 ~ 0*1
i_LS_1 ~ 0*1
i_LS_2 ~ 0*1
i_LS_3 ~ 0*1
i_LS_4 ~ 0*1
i_LS_5 ~ 0*1
i_LS_6 ~ 0*1
i_LS_7 ~ 0*1
t_x ~ 0*1
t_LS ~ 0*1
s_x_1 ~ 0*1
s_x_2 ~ 0*1
s_x_3 ~ 0*1
s_x_4 ~ 0*1
s_x_5 ~ 0*1
s_x_6 ~ 0*1
s_x_7 ~ 0*1
s_LS_1 ~ 0*1
s_LS_2 ~ 0*1
s_LS_3 ~ 0*1
s_LS_4 ~ 0*1
s_LS_5 ~ 0*1
s_LS_6 ~ 0*1
s_LS_7 ~ 0*1
a3 == a2
a4 == a2
a5 == a2
a6 == a2
a7 == a2
b3 == b2
b4 == b2
b5 == b2
b6 == b2
b7 == b2
c3 == c2
c4 == c2
c5 == c2
c6 == c2
c7 == c2
d3 == d2
d4 == d2
d5 == d2
d6 == d2
d7 == d2
xvar2 == xvar1 - a2*a2*xvar1 - d2*d2*yvar1 - 2*a2*cov_ar1*d2
xvar3 == xvar2
xvar4 == xvar2
xvar5 == xvar2
xvar6 == xvar2
xvar7 == xvar2
yvar2 == yvar1 - b2*b2*yvar1 - c2*c2*xvar1 - 2*b2*cov_ar1*c2
yvar3 == yvar2
yvar4 == yvar2
yvar5 == yvar2
yvar6 == yvar2
yvar7 == yvar2
cov_ar2 == (1-a2*b2-c2*d2)*cov_ar1-a2*c2*xvar1-b2*d2*yvar1
cov_ar3 == cov_ar2
cov_ar4 == cov_ar2
cov_ar5 == cov_ar2
cov_ar6 == cov_ar2
cov_ar7 == cov_ar2
sx2 == sx1
sx3 == sx1
sx4 == sx1
sx5 == sx1
sx6 == sx1
sx7 == sx1
sy2 == sy1
sy3 == sy1
sy4 == sy1
sy5 == sy1
sy6 == sy1
sy7 == sy1
xvar1 > 0
yvar1 > 0
x_tVar > 0
y_tVar > 0
cov_txty < (99/100)*sqrt(x_tVar)*sqrt(y_tVar)
cov_txty > -(99/100)*sqrt(x_tVar)*sqrt(y_tVar)
sx1 > 0
sy1 > 0
'

testOut <- lavaan::lavaan(testModel,
                          data2[1:500,],
                          meanstructure = TRUE,
                          missing = 'fiml',
                          int.ov.free = TRUE,
                          int.lv.free = FALSE)


testOut2 <- sem(testModel,
                data2[1:500,])


stab <- .9
xy_r <- .7
cl <- .1

resid <- 1 - stab ^2 - cl^2 - 2*stab*xy_r*cl

xy_r - xy_r * (stab^2) - xy_r * cl ^2 - 2* stab * cl * resid




cor_xyr <- cor_xy - cor_xy * (stab_x * stab_y) -
            cor_xy * (xy * yx) -
            (stab_x * yx * wxr) -
            (stab_y * xy * wyr)

data <- gen_starts(nwaves = 10,
                   ri_x = 1,
                   ri_y = 1,
                   stab_x = .7,
                   stab_y = .7,
                   yx = .00,
                   xy = .00,
                   cor_xy = .5,
                   xr = 0,
                   yr = 0)


testModel <- panelcoder(data2, panelModel="clpm", lags = 3, program = "mplus")
panelPlot(testModel)

testModel5 <- panelcoder(data2, panelModel="clpm", lags = 5, program = "mplus")
panelPlot(testModel5)

testModel9 <- panelcoder(data2, panelModel="clpm", lags = 9, program = "mplus")
panelPlot(testModel9)



waves <- 10
indicators <- 3

data <- gen_starts(
    n = 10000, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .7, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .7, # Stability of X
    stab_y = .7, # Stability of Y
    yx = 0, # Cross lag (Y regressed on X)
    xy = 0, # Cross lag (X regressed on Y)
    cor_xy = .7, # Correlation between X and Y (as correlation)
    xr = 1, # Measurement error for X
    yr = 1 # Measurement error for Y
)

## Create data frame with correctly labeled variables
data2 <- data
names(data2) <- paste(rep(c("x", "LS"), each = waves),
                      rep(1:waves, 2),
                      sep="_")
for (i in names(data2)) {
    data2 <- addIndicators(data2, i, indicators)
}

startC <- (2*waves+1)
endC <- (2*waves)+(2*waves*indicators)
dataI <- data2[,startC:endC]
data2 <- data2[,1:(2*waves)]


## Three Wave

clpm3 <- panelcoder(data2[c(1:3,11:13)],
                    panelModel="clpm",
                    lags = 2,
                    program = "mplus",
                    title = "clpm3")
panelPlot(clpm3)


riclpm3 <- panelcoder(data2[c(1:3,11:13)],
                    panelModel="riclpm",
                    program = "mplus",
                    title = "riclpm3")
panelPlot(riclpm3)


arts3 <- panelcoder(data2[c(1:3,11:13)],
                    panelModel="arts",
                    program = "mplus",
                    title = "arts3")
panelPlot(arts3)



waves <- 10
indicators <- 3
N <- 10000

data <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .7, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .7, # Stability of X
    stab_y = .7, # Stability of Y
    yx = .1, # Cross lag (Y regressed on X)
    xy = .15, # Cross lag (X regressed on Y)
    cor_xy = .7, # Correlation between X and Y (as correlation)
    xr = 1, # Measurement error for X
    yr = 1 # Measurement error for Y
)

## Create data frame with correctly labeled variables
data2 <- data
names(data2) <- paste(rep(c("x", "LS"), each = waves),
                      rep(1:waves, 2),
                      sep="_")
for (i in names(data2)) {
    data2 <- addIndicators(data2, i, indicators)
}

startC <- (2*waves+1)
endC <- (2*waves)+(2*waves*indicators)
dataI <- data2[,startC:endC]
data2 <- data2[,1:(2*waves)]


m1 <- panelcoder(data2,
                 panelModel = "starts",
                 program = "mplus",
                 stationarity = "paths",
                 constrainState = TRUE,
                 analysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001; ITERATIONS=50000")



################################################################################
## STARTS Package
################################################################################

library(STARTS)

data(data.starts01a)
s_data <- data.starts01a[c(2:11)]
names(s_data) <- paste(c(rep("e", 5), rep("n", 5)), c(1:5, 1:5), sep = "_")

test <- panelcoder(s_data[1:5],
                   panelModel = "starts",
                   program = "mplus",
                   stationarity = "full")

s_test <- starts_uni_estimate(s_data[1:5])
summary(s_test)

s_test <- starts_uni_estimate(s_data[1:5], est_var_ar = FALSE)

covmat <- stats::cov( s_data[1:5])
nobs <- nrow(s_data)

mod1a <- starts_uni_estimate(covmat=covmat, nobs=nobs, est_var_ar=FALSE)



test <- panelcoder(data[1:10], panelModel = "dpm_p", slope = "linear", program = "mplus", run = TRUE)


test <- panelcoder(data, panelModel = "alt", slope = "linear", program = "mplus", run = TRUE)


test <- panelcoder(data2[c(1:6,11:16)], panelModel = "gclm", program = "mplus", ma = TRUE, clma = FALSE)

test <- panelcoder(data2[c(1:6)], panelModel = "gclm", program = "mplus", ma = TRUE, clma = FALSE)


test <- panelcoder(data[1:6], panelModel = "gclm", program = "mplus", ma = TRUE, clma = TRUE)


check_nonpos <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
  warning_text <- paste(
    "THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES MAY NOT BE",
    "TRUSTWORTHY FOR SOME PARAMETERS DUE TO A NON-POSITIVE DEFINITE",
    "FIRST-ORDER DERIVATIVE PRODUCT MATRIX.  THIS MAY BE DUE TO THE STARTING",
    "VALUES BUT MAY ALSO BE AN INDICATION OF MODEL NONIDENTIFICATION.",
    sep = " "
  )
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}


check_se <- function(file_path) {
  # Read the entire file as a single string
  text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  
  # Define the target text (collapse newlines and spacing)
  warning_text <- paste(
    "THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES COULD NOT BE",
    "COMPUTED.  THE MODEL MAY NOT BE IDENTIFIED.  CHECK YOUR MODEL.",
    sep = " "
  )
  
  # Normalize whitespace in both strings before searching
  normalize_ws <- function(x) gsub("\\s+", " ", x)
  
  found <- grepl(normalize_ws(warning_text), normalize_ws(text), fixed = TRUE)
    
  return(found)
}




waves <- 10
indicators <- 1
N <- 10000

data_t <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 0, # Random intercept variance for X
    ri_y = 0, # Random intercept variance for Y
    cor_i = .7, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .9, # Stability of X
    stab_y = .9, # Stability of Y
    yx = 0, # Cross lag (Y regressed on X)
    xy = 0, # Cross lag (X regressed on Y)
    cor_xy = .7, # Correlation between X and Y (as correlation)
    xr = 1, # Measurement error for X
    yr = 1 # Measurement error for Y
)

test <- panelcoder(data_t, program = "mplus", limits = TRUE)

check_nonpos("mplus/panelcoder.out")
check_ident("mplus/panelcoder.out")


compareModels(data_t[1:10])

temp <- panelcoder(data2[1:10], program = "mplus", mplusAnalysis = "MODEL=NOCOVARIANCES;\nCOVERAGE=.001;\nITERATIONS=20000;")


data_t <- read_csv("../hildastability/data/ncfsee_w.csv")

test <- panelcoder(data_t[-1], program="lavaan", stationarity = "full")


waves <- 10
indicators <- 1
N <- 1000
cl <- 0
cl2 <- .27

data1 <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 0, # Random intercept variance for X
    ri_y = 0, # Random intercept variance for Y
    cor_i = .7, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = cl, # Cross lag (Y regressed on X)
    xy = cl, # Cross lag (X regressed on Y)
    cor_xy = .3, # Correlation between X and Y (as correlation)
    xr = 0, # Measurement error for X
    yr = 0 # Measurement error for Y
)

data2 <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 0, # Random intercept variance for X
    ri_y = 0, # Random intercept variance for Y
    cor_i = .7, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = cl2, # Cross lag (Y regressed on X)
    xy = cl2, # Cross lag (X regressed on Y)
    cor_xy = .3, # Correlation between X and Y (as correlation)
    xr = 0, # Measurement error for X
    yr = 0 # Measurement error for Y
)


m1 <- panelcoder(data1, program = "mplus", panelModel = "clpm")
m2 <- panelcoder(data2, program = "mplus", panelModel = "clpm")

p.data <- cbind(m1[[5]][1], m2[[5]][1])
names(p.data) <- c("Implied X","Observed X")

plotCors(p.data)



waves <- 10
indicators <- 1
N <- 1000

results <- list()

set.seed(1234)

for (i in 1:20) {
    data <- gen_starts(
        n = N, # N to generate
        nwaves = waves, # Number of waves
        ri_x = 1, # Random intercept variance for X
        ri_y = 1, # Random intercept variance for Y
        cor_i = .5, # Correlation between intercepts (as correlation)
        x = 1, # AR variance for X
        y = 1, # AR variance for Y
        stab_x = .5, # Stability of X
        stab_y = .5, # Stability of Y
        yx = .07, # Cross lag (Y regressed on X)
        xy = .07, # Cross lag (X regressed on Y)
        cor_xy = .5, # Correlation between X and Y (as correlation)
        xr = 2, # Measurement error for X
        yr = 2 # Measurement error for Y
    )
    ## Create data frame with correctly labeled variables
    names(data) <- paste(rep(c("x", "LS"), each = waves),
        rep(1:waves, 2),
        sep = "_"
    )
    out <- panelcoder(data, panelModel = "starts", program = "mplus")
    results[[i]] <- list(
        out[[1]]$xOnY, out[[1]]$xOnY.u,
        out[[1]]$yOnX, out[[1]]$yOnX.u,
        out[[7]], out[[8]]
    )
}
results_df <- map_dfr(
  results,
  ~ tibble(
    v1 = .x[[1]],
    v2 = .x[[2]],
    v3 = .x[[3]],
    v4 = .x[[4]],
    c1 = if (is.null(.x[[5]])) NA_character_ else as.character(.x[[5]]),
    c2 = if (is.null(.x[[6]])) NA_character_ else as.character(.x[[6]])
  )
)
psych::describe(results_df)




flattened <- lapply(results, function(inner_list) {
    lapply(inner_list, function(x) x[[1]])
})

cleaned <- lapply(flattened, function(row) {
    lapply(row, function(x) {
        if (is.null(x)) NA else x
    })
})


results_df <- do.call(rbind.data.frame, cleaned)
    

set.seed(123)

data <- gen_starts(
    n = N, # N to generate
    nwaves = waves, # Number of waves
    ri_x = 1, # Random intercept variance for X
    ri_y = 1, # Random intercept variance for Y
    cor_i = .5, # Correlation between intercepts (as correlation)
    x = 1, # AR variance for X
    y = 1, # AR variance for Y
    stab_x = .5, # Stability of X
    stab_y = .5, # Stability of Y
    yx = .07, # Cross lag (Y regressed on X)
    xy = .07, # Cross lag (X regressed on Y)
    cor_xy = .7, # Correlation between X and Y (as correlation)
    xr = 2, # Measurement error for X
    yr = 2 # Measurement error for Y
)

    
## Create data frame with correctly labeled variables
names(data) <- paste(rep(c("x", "LS"), each = waves),
    rep(1:waves, 2),
    sep = "_"
)
out <- panelcoder(data, panelModel = "starts", program = "mplus")

out <- panelcoder(data, panelModel = "riclpm", program = "mplus")

out <- panelcoder(data, panelModel = "arts", program = "mplus")



data <- gen_starts(n=10000,
                   nwaves = 5,
                   ri_x = 1,
                   ri_y = 1,
                   stab_x = .5,
                   stab_y = .5,
                   yx = .07,
                   xy = .07,
                   cor_i = .33,
                   cor_xy = .33,
                   xr = 0,
                   yr = 0)

mr <- panelcoder(data, panelModel="riclpm", program = "mplus")
md <- panelcoder(data, panelModel="dpm_p", program = "mplus")


waves <- 5
lagss <- 4

test <- panelcoder(data2m[c(1:5,11:15)],
                   panelModel="clpm",
                   program = "lavaan",
                   stationarity = "paths",
                   lags = 1)


data2m <- data2
data2m[c(1,3,5,9,11),] <- NA
