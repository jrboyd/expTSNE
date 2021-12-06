testthat::context("expTSNE creation and validity")
library(expTSNE)
library(testthat)

test_counts = matrix(runif(10000), nrow = 10, ncol = 1000)
rownames(test_counts) = paste("row", seq_len(nrow(test_counts)))
colnames(test_counts) = paste("col", seq_len(ncol(test_counts)))

test_that("expTSNE.input minimal", {
  et = expTSNE.input(test_counts)
  testthat::expect_is(et, "expTSNE.input")
})

test_that("expTSNE.input need rownames", {
  bad_counts = test_counts
  rownames(bad_counts) = NULL
  testthat::expect_error({et = expTSNE.input(bad_counts)}, regexp = "raw_counts must have rownames.")
})

test_that("expTSNE.input need colnames", {
  bad_counts = test_counts
  colnames(bad_counts) = NULL
  testthat::expect_error({et = expTSNE.input(bad_counts)}, regexp = "raw_counts must have colnames.")
})

test_that("expTSNE.input selected_rows in rownames", {
  testthat::expect_error({et = expTSNE.input(test_counts, selected_rows = "not")}, regexp = "selected_rows must all be in rownames.")
})

test_that("expTSNE.input selected_columns in colnames", {
  testthat::expect_error({et = expTSNE.input(test_counts, selected_columns = "not")}, regexp = "selected_columns must all be in colnames.")
})

test_that("expTSNE.input raw_counts and norm_counts match", {
  testthat::expect_error({et = expTSNE.input(test_counts, norm_counts = test_counts[-1,])}, regexp = "norm_counts and raw_counts rownames must be equal.")
  testthat::expect_error({et = expTSNE.input(test_counts, norm_counts = test_counts[,-1])}, regexp = "norm_counts and raw_counts colnames must be equal.")
})

test_that("expTSNE.input metadata must have column_id", {
  # et = testthat::expect_error({expTSNE.input(test_counts, meta_data = data.frame(row.names = NULL))}, regexp = "meta_data must have rownames.")
  testthat::expect_error({expTSNE.input(test_counts, meta_data = data.frame(1))}, regexp = "meta_data must have one row per column in raw_counts.")
  mdf = data.frame(seq(ncol(test_counts)))
  testthat::expect_error({expTSNE.input(test_counts, meta_data = mdf)}, regexp = "meta_data rownames must equal raw_counts colnames.")
  rownames(mdf) = colnames(test_counts)
  et = expTSNE.input(test_counts, meta_data = mdf)
  testthat::expect_is(et, "expTSNE.input")
})
