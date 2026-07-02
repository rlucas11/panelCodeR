# Layer 2: panelcoder() input-validation guard rails.
#
# These are the stop()/warning() conditions that protect users from asking for
# unsupported model specifications. They run through panelcoder(run = FALSE), so
# they're fast, and they lock in the safety behaviour so future refactors can't
# silently drop a check. Each uses a well-formed fixture so failure is due to
# the specification, not the data.

biv <- make_pc_data(nwaves = 5, indicators = 1, bivariate = TRUE)
uni <- make_pc_data(nwaves = 5, indicators = 1, bivariate = FALSE)

# --- errors -------------------------------------------------------------------

test_that("unknown model name is rejected", {
    expect_error(
        pc_code(biv, panelModel = "nonsense"),
        "No model with that name"
    )
})

test_that("measurement model needs more than one indicator", {
    expect_error(
        pc_code(uni, panelModel = "measurement"),
        "only 1 indicator"
    )
})

test_that("invalid slope option is rejected", {
    expect_error(
        pc_code(biv, panelModel = "starts", slope = "quadratic"),
        "Use 'linear', 'basis', or 'centered'"
    )
})

test_that("lgcm requires an explicit slope", {
    expect_error(
        pc_code(biv, panelModel = "lgcm"),
        "Need to specify"
    )
})

test_that("constrained DPM does not allow a slope", {
    expect_error(
        pc_code(biv, panelModel = "dpm_c", slope = "linear"),
        "Constrained DPM"
    )
})

test_that("lags > 1 not allowed with a random-intercept model", {
    expect_error(
        pc_code(biv, panelModel = "starts", lags = 2),
        "Lags greater than 1"
    )
})

test_that("phantom (missing) waves require full stationarity", {
    phantom <- make_pc_data(waves = c(1, 2, 4))
    expect_error(
        pc_code(phantom, panelModel = "clpm"),
        "phantom variables when full stationarity"
    )
})

# --- warnings -----------------------------------------------------------------
# pc_code() suppresses warnings to keep code snapshots clean, so warning
# assertions call panelcoder() directly and swallow the printed model code.

expect_pc_warning <- function(regexp, ...) {
    expect_warning(
        suppressMessages(utils::capture.output(
            panelcoder(..., run = FALSE)
        )),
        regexp
    )
}

test_that("ALT forces slope to linear with a warning", {
    expect_pc_warning("Slope set to linear", biv, panelModel = "alt")
})

test_that("lags > 1 with full stationarity warns and downgrades", {
    expect_pc_warning(
        "Stationarity not implemented with lags",
        biv, panelModel = "clpm", lags = 2, stationarity = "full"
    )
})
