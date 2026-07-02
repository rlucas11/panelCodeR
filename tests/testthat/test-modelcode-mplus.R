# Layer 3: generated Mplus code (snapshot regression tests).
#
# Mirror of test-modelcode-lavaan.R for program = "mplus". panelcoder(run = FALSE)
# returns the MODEL syntax from lav2mplus() without launching Mplus, so these are
# fast and deterministic. See that file's header for the accept/review workflow.

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
        test_that(paste0("mplus code: ", model, " (bivariate, single indicator)"), {
            args <- c(list(data = biv, panelModel = model, program = "mplus"), extra)
            expect_snapshot(cat(do.call(pc_code, args)))
        })
    })
}

test_that("mplus code: measurement model (bivariate, 3 indicators)", {
    d <- make_pc_data(nwaves = 4, indicators = 3, bivariate = TRUE)
    expect_snapshot(cat(pc_code(d, panelModel = "measurement", program = "mplus")))
})

test_that("mplus code: starts (univariate, single indicator)", {
    d <- make_pc_data(nwaves = 5, indicators = 1, bivariate = FALSE)
    expect_snapshot(cat(pc_code(d, panelModel = "starts", program = "mplus")))
})

test_that("mplus code: starts (bivariate, 3 indicators)", {
    d <- make_pc_data(nwaves = 5, indicators = 3, bivariate = TRUE)
    expect_snapshot(cat(pc_code(d, panelModel = "starts", program = "mplus")))
})
