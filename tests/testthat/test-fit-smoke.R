# Layer 6: end-to-end smoke tests (SLOW -- these actually fit models).
#
# Everything else in the suite is fast because nothing is estimated. These few
# tests round it out by confirming that the generated lavaan code is *accepted*
# by lavaan and that panelcoder() returns a well-formed pcOutput. They are gated
# with skip_on_cran() so the normal, fast `devtools::test()` on CRAN/CI stays
# quick; run them locally when you want end-to-end coverage.
#
# They deliberately assert only the return contract and that fitting produced no
# fatal error -- estimate *correctness* is lavaan's job, not panelCodeR's, so we
# don't test parameter values (which would be slow and brittle).

# STARTS-structured data so the models are actually estimable. Multi-indicator
# variant (indicators = 3) is used for the measurement model.
make_starts_data <- function(nwaves = 6, n = 600, indicators = 1,
                            seed = 20240702) {
    set.seed(seed)
    d <- gen_starts(n = n, nwaves = nwaves,
                    ri_x = 1, ri_y = 1, cor_i = .5,
                    x = 1, y = 1, stab_x = .5, stab_y = .5,
                    yx = .2, xy = .1, cor_xy = .4, xr = 1, yr = 1)
    names(d) <- paste(rep(c("x", "y"), each = nwaves),
                      rep(seq_len(nwaves), 2), sep = "_")
    if (indicators > 1) {
        for (v in names(d)) d <- addIndicators(d, v, indicators)
        d <- d[, (2 * nwaves + 1):ncol(d)]   # keep only the indicator columns
    }
    d
}

# Fit (suppressing the printed summary) and assert the pcOutput contract.
expect_fits <- function(data, ...) {
    out <- NULL
    suppressWarnings(utils::capture.output(
        out <- panelcoder(data = data, program = "lavaan", run = TRUE, ...)
    ))
    expect_s3_class(out, "pcOutput")
    expect_length(out, 8)
    expect_false(is.null(out[[1]]), label = "pcSum (summary) present on success")
    expect_null(out[[8]], label = "errorM (fatal error)")
    invisible(out)
}

smoke_models <- list(
    clpm   = list(),
    riclpm = list(),
    starts = list(),
    arts   = list(),
    sts    = list(),
    dpm_c  = list(),
    dpm_p  = list(),
    gclm   = list(),
    lgcm   = list(slope = "linear")
)

biv <- make_starts_data()

for (m in names(smoke_models)) {
    local({
        model <- m
        extra <- smoke_models[[m]]
        test_that(paste0(model, " fits end to end (lavaan)"), {
            skip_on_cran()
            skip_if_not_installed("lavaan")
            do.call(expect_fits, c(list(data = biv, panelModel = model), extra))
        })
    })
}

test_that("measurement model fits end to end (lavaan)", {
    skip_on_cran()
    skip_if_not_installed("lavaan")
    d <- make_starts_data(indicators = 3)
    expect_fits(d, panelModel = "measurement")
})
