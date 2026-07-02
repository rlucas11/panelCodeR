# Layer 1: getInfo() name parsing and validation.
#
# getInfo() is the foundation of the whole pipeline: it turns column names into
# the `info` structure every downstream function relies on. It's pure and fast,
# so we test both correct parsing and its guard rails here.

test_that("univariate single-indicator names parse correctly", {
    info <- getInfo(make_pc_data(nwaves = 5, indicators = 1, bivariate = FALSE))

    expect_false(info$gen$yVar)
    expect_equal(info$gen$maxWaves, 5)
    expect_equal(info$x$name, "x")
    expect_equal(info$x$indicators, 1)
    expect_equal(info$x$waves, 1:5)
    expect_equal(info$x$actualWaves, 1:5)
    expect_null(info$y)
})

test_that("univariate multi-indicator names parse correctly", {
    info <- getInfo(make_pc_data(nwaves = 4, indicators = 3, bivariate = FALSE))

    expect_false(info$gen$yVar)
    expect_equal(info$gen$maxWaves, 4)
    expect_equal(info$x$indicators, 3)
})

test_that("bivariate names set yVar and the second variable", {
    info <- getInfo(make_pc_data(nwaves = 5, indicators = 1, bivariate = TRUE))

    expect_true(info$gen$yVar)
    expect_equal(info$x$name, "x")
    expect_equal(info$y$name, "y")
    expect_equal(info$y$indicators, 1)
})

test_that("bivariate multi-indicator names parse correctly", {
    info <- getInfo(make_pc_data(nwaves = 4, indicators = 3, bivariate = TRUE))

    expect_true(info$gen$yVar)
    expect_equal(info$x$indicators, 3)
    expect_equal(info$y$indicators, 3)
})

test_that("missing waves are detected as phantom (actualWaves != waves)", {
    info <- getInfo(make_pc_data(waves = c(1, 2, 4)))

    expect_equal(info$gen$maxWaves, 4)
    expect_equal(info$x$waves, 1:4)          # full, contiguous range
    expect_equal(info$x$actualWaves, c(1, 2, 4))  # what's actually present
})

# --- validation / guard rails -------------------------------------------------

named_df <- function(nms) {
    df <- as.data.frame(matrix(0, nrow = 2, ncol = length(nms)))
    names(df) <- nms
    df
}

test_that("more than two variable stems is an error", {
    expect_error(
        getInfo(named_df(c("x_1", "y_1", "z_1"))),
        "more than two variables"
    )
})

test_that("too many name splits is an error", {
    expect_error(
        getInfo(named_df(c("x_1_1_1", "x_1_1_2"))),
        "Too may splits"
    )
})

test_that("inconsistent splitting within a variable is an error", {
    expect_error(
        getInfo(named_df(c("x_1", "x_1_1"))),
        "cannot be split in the same way"
    )
})

test_that("differing indicator counts across waves is an error", {
    expect_error(
        getInfo(named_df(c("x_1_1", "x_1_2", "x_2_1"))),
        "Different number of indicators per wave"
    )
})

test_that("names longer than 8 characters warn (Mplus limit)", {
    expect_warning(
        getInfo(named_df(c("verylongname_1", "verylongname_2"))),
        "longer than 8 characters"
    )
})
