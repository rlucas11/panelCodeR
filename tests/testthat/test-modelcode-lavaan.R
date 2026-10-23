# Layer 3: generated lavaan code (snapshot regression tests).
#
# These snapshot the exact model syntax panelcoder() produces for each model,
# without ever fitting anything. When you change code generation, `devtools::test()`
# will show a readable diff; review it and run `testthat::snapshot_accept()` to
# bless intended changes. Whether lavaan then *fits* the code is out of scope --
# that's lavaan's job, not panelCodeR's.
#
# The main matrix is bivariate, single-indicator, 5 waves. Model families with
# special requirements carry their needed argument in `model_specs`. Measurement
# and multi-indicator variants get dedicated tests below.

model_specs <- list(
    starts = list(),
    riclpm = list(),
    arts   = list(),
    clpm   = list(),
    sts    = list(),
    dpm_c  = list(),
    dpm_p  = list(),
    gclm   = list(),
    lgcm   = list(slope = "linear"),
    alt    = list(slope = "linear"),
    lcmsr  = list(slope = "linear")
)

biv <- make_pc_data(nwaves = 5, indicators = 1, bivariate = TRUE)

for (m in names(model_specs)) {
    local({
        model <- m
        extra <- model_specs[[m]]
        test_that(paste0("lavaan code: ", model, " (bivariate, single indicator)"), {
            args <- c(list(data = biv, panelModel = model, program = "lavaan"), extra)
            expect_snapshot(cat(do.call(pc_code, args)))
        })
    })
}

test_that("lavaan code: measurement model (bivariate, 3 indicators)", {
    d <- make_pc_data(nwaves = 4, indicators = 3, bivariate = TRUE)
    expect_snapshot(cat(pc_code(d, panelModel = "measurement", program = "lavaan")))
})

test_that("lavaan code: starts (univariate, single indicator)", {
    d <- make_pc_data(nwaves = 5, indicators = 1, bivariate = FALSE)
    expect_snapshot(cat(pc_code(d, panelModel = "starts", program = "lavaan")))
})

test_that("lavaan code: starts (bivariate, 3 indicators)", {
    d <- make_pc_data(nwaves = 5, indicators = 3, bivariate = TRUE)
    expect_snapshot(cat(pc_code(d, panelModel = "starts", program = "lavaan")))
})

test_that("lavaan code: starts with free trait loadings (bivariate, single indicator)", {
    expect_snapshot(cat(pc_code(biv, panelModel = "starts", traitLoadings = "free",
                                program = "lavaan")))
})

test_that("free trait loadings stay fixed on phantom waves", {
    phantom <- make_pc_data(waves = c(1, 2, 4, 5))
    model <- suppressWarnings(.buildModel(phantom, panelModel = "starts",
                                          stationarity = "full",
                                          traitLoadings = "free"))
    loadings <- model[model$lhs == "t_x" & model$op == "=~", ]
    expect_equal(loadings$free[loadings$rhs %in% c("l_x_1", "l_x_3")], c(0, 0))
    expect_equal(loadings$ustart[loadings$rhs %in% c("l_x_1", "l_x_3")], c(1, 1))
    expect_true(all(loadings$free[loadings$rhs %in% c("l_x_2", "l_x_4", "l_x_5")] == 1))
})
