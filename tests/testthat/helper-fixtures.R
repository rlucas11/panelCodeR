# Helpers shared across test files.
#
# Because panelCodeR builds its parameter table (and therefore all generated
# model code) purely from the *column names* of the data -- getInfo() parses
# names, and the build path never touches data values -- test fixtures only
# need correctly named columns plus a few dummy rows. This keeps the code-
# generation tests fully deterministic and independent of gen_starts().

# Build a minimal, correctly-named data frame for a given design.
#
#   nwaves     number of waves
#   indicators indicators per wave (1 = single-indicator)
#   bivariate  add a second variable stem
#   waves      optional explicit wave indices (e.g. c(1, 2, 4) to create a
#              phantom / missing wave); defaults to 1:nwaves
#   stems      variable stems to use
make_pc_data <- function(nwaves = 5,
                         indicators = 1,
                         bivariate = FALSE,
                         waves = NULL,
                         stems = c("x", "y"),
                         nobs = 25,
                         seed = 1) {
    if (is.null(waves)) waves <- seq_len(nwaves)

    build_stem <- function(stem) {
        if (indicators == 1) {
            paste(stem, waves, sep = "_")
        } else {
            # wave-major, indicators nested within wave: x_1_1, x_1_2, x_2_1 ...
            as.vector(t(outer(waves, seq_len(indicators),
                              function(w, i) paste(stem, w, i, sep = "_"))))
        }
    }

    vars <- build_stem(stems[1])
    if (bivariate) vars <- c(vars, build_stem(stems[2]))

    set.seed(seed)
    df <- as.data.frame(matrix(stats::rnorm(nobs * length(vars)), nrow = nobs))
    names(df) <- vars
    df
}

# Generate model code without fitting. panelcoder(run = FALSE) cat()s the code
# and returns it; we discard the printed side-effect and hand back the string.
# Expected warnings (e.g. auto-setting a slope) are suppressed so code snapshots
# stay focused on the code itself -- warnings are asserted separately in
# test-validation.R.
pc_code <- function(data, ...) {
    out <- NULL
    suppressWarnings(
        utils::capture.output(out <- panelcoder(data = data, ..., run = FALSE))
    )
    out
}
