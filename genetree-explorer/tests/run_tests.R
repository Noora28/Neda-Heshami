#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  tests_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  tests_dir <- normalizePath("tests")
}
app_dir <- dirname(tests_dir)
setwd(app_dir)

source(file.path(app_dir, "R", "ncbi_api.R"), local = FALSE)
source(file.path(app_dir, "R", "phylo.R"), local = FALSE)
.test_dir <<- tests_dir

reporter <- testthat::ProgressReporter$new()
results <- testthat::test_dir(tests_dir, reporter = reporter, wrap = TRUE)
df <- as.data.frame(results)
n_fail <- if ("failed" %in% names(df)) sum(df$failed) else 0
n_err <- if ("error" %in% names(df)) sum(df$error) else 0
if (n_fail + n_err > 0) {
  message("Fixture tests failed.")
  quit(status = 1)
}
message("All fixture tests passed.")
