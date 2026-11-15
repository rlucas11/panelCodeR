# Starting values (startVals argument). These tests only build parameter
# tables and model code; nothing is fitted.

biv <- make_pc_data(nwaves = 5, indicators = 1, bivariate = TRUE, nobs = 200)

build <- function(data, startVals, ...) {
    model <- suppressWarnings(.buildModel(data, ...))
    .applyStartVals(model, .resolveStartVals(startVals, model, getInfo(data), data))
}

ustart_of <- function(model, label) model$ustart[model$label == label]

test_that("startVals = NULL leaves the table unchanged", {
    model <- suppressWarnings(.buildModel(biv, panelModel = "starts"))
    expect_identical(build(biv, NULL, panelModel = "starts"), model)
})

test_that("starts are only set on free, non-dependent parameters", {
    model <- build(biv, "data", panelModel = "starts", stationarity = "full")
    deps <- .dependentLabels(model)
    set <- !is.na(model$ustart) & model$free != 0
    expect_false(any(model$label[set] %in% deps))
    # fixed parameters keep their fixed values
    base <- suppressWarnings(.buildModel(biv, panelModel = "starts",
                                         stationarity = "full"))
    expect_identical(model$ustart[model$free == 0], base$ustart[base$free == 0])
    # dependent parameters under full stationarity stay unset
    expect_true(is.na(ustart_of(model, "xvar2")))
    expect_true(is.na(ustart_of(model, "cov_ar2")))
    expect_true(is.na(ustart_of(model, "a_1_3")))
})

test_that("data rule splits variance across components and uses .5 / 0", {
    model <- build(biv, "data", panelModel = "starts", stationarity = "full")
    v <- mean(sapply(paste0("x_", 1:5), function(n) var(biv[[n]])))
    expect_equal(ustart_of(model, "x_tVar"), signif(v / 3, 4))
    expect_equal(ustart_of(model, "xvar1"), signif(v / 3, 4))
    expect_equal(ustart_of(model, "sx1"), signif(v / 3, 4))
    expect_equal(ustart_of(model, "a_1_2"), .5)
    expect_equal(ustart_of(model, "c_1_2"), 0)
    expect_equal(ustart_of(model, "cov_txty"), 0)

    riclpm <- build(biv, "data", panelModel = "riclpm", stationarity = "full")
    expect_equal(ustart_of(riclpm, "x_tVar"), signif(v / 2, 4))
    clpm <- build(biv, "data", panelModel = "clpm", stationarity = "full")
    expect_equal(ustart_of(clpm, "xvar1"), signif(v, 4))
})

test_that("data rule gives innovation variances when not stationary", {
    model <- build(biv, "data", panelModel = "riclpm", stationarity = "paths")
    v2 <- var(biv$x_2)
    expect_equal(ustart_of(model, "xvar1"), signif(var(biv$x_1) / 2, 4))
    expect_equal(ustart_of(model, "xvar2"), signif(v2 / 2 * .75, 4))
})

test_that("data rule covers every free parameter across models", {
    specs <- list(
        list(panelModel = "starts"), list(panelModel = "riclpm"),
        list(panelModel = "arts"), list(panelModel = "clpm"),
        list(panelModel = "sts"), list(panelModel = "dpm_c"),
        list(panelModel = "dpm_p"),
        list(panelModel = "gclm", ma = TRUE, clma = TRUE),
        list(panelModel = "lgcm", slope = "basis"),
        list(panelModel = "alt", slope = "linear"),
        list(panelModel = "lcmsr", slope = "linear"),
        list(panelModel = "clpm", lags = 2),
        list(panelModel = "starts", traitLoadings = "free")
    )
    for (ind in c(1, 3)) {
        d <- make_pc_data(nwaves = 4, indicators = ind, bivariate = TRUE,
                          nobs = 200)
        for (s in specs) for (st in c("none", "paths", "full")) {
            model <- suppressWarnings(do.call(.buildModel,
                                              c(list(data = d, stationarity = st), s)))
            vals <- .resolveStartVals("data", model, getInfo(d), d)
            el <- .startValEligible(model)
            expect_true(all(is.finite(vals[el])),
                        info = paste(s$panelModel, st, ind))
        }
    }
})

test_that("labels override specific roles, which override general roles", {
    model <- build(biv, list(stab = .7, stab_y = .6, a_1_2 = .3, cl_xy = .1),
                   panelModel = "riclpm", stationarity = "none")
    expect_equal(ustart_of(model, "a_1_2"), .3)   # label
    expect_equal(ustart_of(model, "a_1_3"), .7)   # general role
    expect_equal(ustart_of(model, "b_1_2"), .6)   # specific role
    expect_equal(ustart_of(model, "c_1_2"), .1)   # x -> y
    expect_true(is.na(ustart_of(model, "d_1_2"))) # y -> x not given
})

test_that("unknown and constrained names warn", {
    model <- suppressWarnings(.buildModel(biv, panelModel = "starts",
                                          stationarity = "full"))
    info <- getInfo(biv)
    expect_warning(.resolveStartVals(list(bogus = 1), model, info, biv),
                   "not roles or parameter labels")
    expect_warning(.resolveStartVals(c(xvar2 = 1), model, info, biv),
                   "determined by constraints")
    expect_error(.resolveStartVals(list(.5), model, info, biv))
    expect_error(.resolveStartVals("defaults", model, info, biv))
})

test_that("a function can supply start values", {
    model <- build(biv, function(data, info) list(stab = .4),
                   panelModel = "clpm")
    expect_equal(ustart_of(model, "a_1_2"), .4)
})

test_that("generated code contains start values", {
    mplus <- pc_code(biv, panelModel = "starts", program = "mplus",
                     stationarity = "full", startVals = list(stab = .5, cl = 0))
    expect_match(mplus, "a_x_2 ON a_x_1*0.5 (a_1_2);", fixed = TRUE)
    expect_match(mplus, "a_x_2 ON a_y_1*0 (d_1_2);", fixed = TRUE)
    expect_match(mplus, "a_x_3 ON a_x_2* (a_1_3);", fixed = TRUE)

    lav <- pc_code(biv, panelModel = "starts", program = "lavaan",
                   stationarity = "full", startVals = list(stab = .5))
    expect_match(lav, "a_x_2 ~ start(0.5)*a_x_1 + a_1_2*a_x_1", fixed = TRUE)
    expect_match(lav, "a_x_3 ~ a_1_3*a_x_2", fixed = TRUE)
})

test_that("small and large values are written without scientific notation", {
    expect_equal(.formatStartVal(c(.00001234, 123456, .5)),
                 c("0.00001234", "123456", "0.5"))
})

test_that("start value sequences are recognised", {
    expect_true(.isStartValsSeq(list(NULL, "data")))
    expect_false(.isStartValsSeq(list(stab = .5)))
    expect_false(.isStartValsSeq("data"))
    expect_error(pc_code(biv, startVals = list(NULL, "data")), "run = TRUE")
    expect_error(pc_code(biv, startValsRule = "last"), "startValsRule")
})
