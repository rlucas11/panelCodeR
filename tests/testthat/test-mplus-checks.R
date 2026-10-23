# check_mplus_error() reads an Mplus .out file and reports whether Mplus failed
# before estimation, so panelcoder() can return a clear error message instead of
# failing later when summarizing missing results.

write_out <- function(lines) {
    path <- tempfile(fileext = ".out")
    writeLines(lines, path)
    path
}

test_that("missing output file is reported", {
    expect_match(check_mplus_error(tempfile(fileext = ".out")),
                 "output file not found")
})

test_that("Mplus input errors are reported with their explanation", {
    path <- write_out(c(
        "INPUT READING TERMINATED NORMALLY",
        "",
        "*** ERROR in DATA command",
        "  The file specified for the FILE option cannot be found.  Check that this",
        "  file exists:  mydata.dat",
        "",
        "MUTHEN & MUTHEN"
    ))
    msg <- check_mplus_error(path)
    expect_match(msg, "^\\*\\*\\* ERROR in DATA command")
    expect_match(msg, "file exists:\\s+mydata.dat$")
    expect_no_match(msg, "MUTHEN")
})

test_that("output without errors returns NULL", {
    path <- write_out(c("INPUT READING TERMINATED NORMALLY", "", "MODEL RESULTS"))
    expect_null(check_mplus_error(path))
})
