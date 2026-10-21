# predetermined = TRUE should only correlate initial impulses with factors
# (stable trait, slope) that actually exist in the model.

biv <- make_pc_data(nwaves = 5, indicators = 1, bivariate = TRUE)

predet_lines <- function(code) {
    lines <- trimws(strsplit(code, "\n")[[1]])
    lines[grepl("^(t|sl)_\\w+ ~~ .*i_\\w+_1$", lines)]
}

test_that("no trait covariances in models without a stable trait", {
    for (m in c("arts", "clpm")) {
        code <- pc_code(biv, panelModel = m, predetermined = TRUE)
        expect_length(predet_lines(code), 0)

        code <- pc_code(biv, panelModel = m, predetermined = TRUE,
                        slope = "linear")
        lines <- predet_lines(code)
        expect_length(lines, 4)
        expect_false(any(startsWith(lines, "t_")))
    }
})

test_that("models with a stable trait keep trait and slope covariances", {
    code <- pc_code(biv, panelModel = "riclpm", predetermined = TRUE,
                    slope = "linear")
    lines <- predet_lines(code)
    expect_length(lines, 8)
    expect_true("sl_y ~~ c_sly_ix*i_x_1" %in% lines)
    expect_true("sl_y ~~ c_sly_iy*i_y_1" %in% lines)
})

test_that("predetermined is ignored in models without impulses", {
    code <- pc_code(biv, panelModel = "sts", predetermined = TRUE)
    expect_length(predet_lines(code), 0)
})

test_that("dpm_p with full stationarity leaves wave 1 free", {
    expect_warning(
        code <- panelcoder(biv, panelModel = "dpm_p", stationarity = "full",
                           run = FALSE) |>
            utils::capture.output() |> invisible(),
        "incompatible with the predetermined DPM"
    )
    code <- pc_code(biv, panelModel = "dpm_p", stationarity = "full")
    lines <- trimws(strsplit(code, "\n")[[1]])
    expect_false(any(grepl("^(xvar2|yvar2|cov_ar2) ==", lines)))
    expect_true("xvar3 == xvar2" %in% lines)
    expect_true("yvar3 == yvar2" %in% lines)
    expect_true("cov_ar3 == cov_ar2" %in% lines)
})
