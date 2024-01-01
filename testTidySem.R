library(panelCodeR)
library(tidyverse)
library(tidySEM)
library(lavaanExtra)


data <- gen_starts(
    n = 10000, # N to generate
    nwaves = 5, # Number of waves
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

################################################################################
## panelCodeR
################################################################################

data2 <- data
names(data) <- c(paste0("x", 1:5), paste0("y", 1:5))

pcr <- run_starts_lavaan(data, 5)

pcr <- lavaanStarts2(5, xWaves=1:5, yWaves=1:5)

lavaanify(pcr)

summary(pcr)

################################################################################
## tidySEM
################################################################################



names(data) <- paste(names(data), 1, sep="_")

model <- tidy_sem(data)

model |>
    measurement(
        auto.var = TRUE,
        auto.cov.lv.x = FALSE,
        auto.efa = FALSE,
        auto.cov.y. = FALSE,
        meanstructure = FALSE
    ) -> model
model

model |>
    add_paths(
        c(
            "ri_x = ~ 1 * x1 + 1 * x2 + 1 * x3 + 1 * x4 + 1 * x5",
            "ri_y = ~ 1 * y1 + 1 * y2 + 1 * y3 + 1 * y4 + 1 * y5",
            "arx1 = ~ 1 * x1",
            "arx2 = ~ 1 * x2",
            "arx3 = ~ 1 * x3",
            "arx4 = ~ 1 * x4",
            "arx5 = ~ 1 * x5",
            "ary1 = ~ 1 * y1",
            "ary2 = ~ 1 * y2",
            "ary3 = ~ 1 * y3",
            "ary4 = ~ 1 * y4",
            "ary5 = ~ 1 * y5",
            "ri_x ~~ x.t.var*ri_x",
            "ri_y ~~ y.t.var*ri_y",
            "arx2 ~ a*arx1 + d*ary1",
            "arx3 ~ a*arx2 + d*ary2",
            "arx4 ~ a*arx3 + d*ary3",
            "arx5 ~ a*arx4 + d*ary4",
            "ary2 ~ b*ary1 + c*arx1",
            "ary3 ~ b*ary2 + c*arx2",
            "ary4 ~ b*ary3 + c*arx3",
            "ary5 ~ b*ary4 + c*arx4",
            "arx1 ~~ x.ar.var*arx1",
            "arx2 ~~ x.ar.var2*arx2",
            "arx3 ~~ x.ar.var2*arx3",
            "arx4 ~~ x.ar.var2*arx4",
            "arx5 ~~ x.ar.var2*arx5",
            "ary1 ~~ y.ar.var*ary1",
            "ary2 ~~ y.ar.var2*ary2",
            "ary3 ~~ y.ar.var2*ary3",
            "ary4 ~~ y.ar.var2*ary4",
            "ary5 ~~ y.ar.var2*ary5",
            "ri_x ~~ ri_y",
            "arx1 ~~ ary1",
            "x1 ~~ 0*x1",
            "x2 ~~ 0*x2",
            "x3 ~~ 0*x3",
            "x4 ~~ 0*x4",
            "x5 ~~ 0*x5",
            "y1 ~~ 0*y1",
            "y2 ~~ 0*y2",
            "y3 ~~ 0*y3",
            "y4 ~~ 0*y4",
            "y5 ~~ 0*y5"
        ),
        meanstructure = FALSE,
        auto.cov.lv.x = FALSE,
        auto.efa = FALSE,
        auto.cov.y = FALSE
    ) -> model

dictionary(model) <- rbind(
    dictionary(model),
    c("ri_x", NA, "latent", "ri_x"),
    c("ri_y", NA, "latent", "ri_y")
)


syntax(model) <- syntax(model)[1:43, ]


test <- model |>
    estimate_lavaan(func = "lavaan", meanstructure = FALSE)


################################################################################
## lavaanExtra
################################################################################

dfNames <- names(data)

latent <- list(
    x1 = dfNames[1])


setNames(as.list(c(1,2)), c("foo", "bar"))
    


################################################################################
## tidySEM Example
################################################################################

df <- HolzingerSwineford1939
names(df)

names(df)[grepl("^x", names(df))] <- c("vis_1", "vis_2", "vis_3", "tex_1", "tex_2", "tex_3", "spe_1", "spe_2", "spe_3")

model <- tidy_sem(df)

model |>
  measurement() -> model

model |>
    add_paths(c(
        "vis =~ 1*vis_1",
        "vis =~ a*vis_2",
        "vis =~ a*vis_3"
    ),
    auto.fix.first = FALSE) -> model

fit <- estimate_lavaan(model)
summary(fit)

model$syntax[35, "free"] <- 1
model$syntax[36, "free"] <- 1


fit2 <- estimate_lavaan(model)
summary(fit)


model$syntax[2, "label"] <- "a"
model$syntax[3, "label"] <- "a"


model$syntax[22, "label"] <- "a"
model$syntax[23, "label"] <- "a"
model$syntax[24, "label"] <- "a"

fit2 <- estimate_lavaan(model)
summary(fit2)

