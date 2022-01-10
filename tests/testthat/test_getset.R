testthat::context("getters and setters")
library(expTSNE)
library(testthat)

test_counts = matrix(runif(10000), nrow = 10, ncol = 1000)
rownames(test_counts) = paste("row", seq_len(nrow(test_counts)))
colnames(test_counts) = paste("col", seq_len(ncol(test_counts)))

pkg_dir = system.file(package = "expTSNE", "extdata")

et.input.loaded = expTSNE.load(file.path(pkg_dir, "test_expTSNE.input"))
et.tsne.loaded = expTSNE.load(file.path(pkg_dir, "test_expTSNE"))

test_that("names", {
  expect_equal(names(et.input.loaded), c("raw_counts", "norm_counts", "meta_data", "perplexity", "seed", "selected_rows", "selected_columns"))
  expect_equal(names(et.tsne.loaded), c("raw_counts", "norm_counts", "meta_data", "perplexity", "seed", "selected_rows", "selected_columns", 'tsne_result'))
})

test_that("get input", {
  expect_equal(et.input.loaded$raw_counts, et.input.loaded@raw_counts)
  expect_equal(et.input.loaded$norm_counts, et.input.loaded@norm_counts)
  expect_equal(et.input.loaded$meta_data, et.input.loaded@meta_data)
  expect_equal(et.input.loaded$perplexity, et.input.loaded@perplexity)
  expect_equal(et.input.loaded$seed, et.input.loaded@seed)
  expect_equal(et.input.loaded$selected_rows, et.input.loaded@selected_rows)
  expect_equal(et.input.loaded$selected_columns, et.input.loaded@selected_columns)
})

test_that("get tsne", {
  expect_equal(et.tsne.loaded$raw_counts, et.tsne.loaded@raw_counts)
  expect_equal(et.tsne.loaded$norm_counts, et.tsne.loaded@norm_counts)
  expect_equal(et.tsne.loaded$meta_data, et.tsne.loaded@meta_data)
  expect_equal(et.tsne.loaded$perplexity, et.tsne.loaded@perplexity)
  expect_equal(et.tsne.loaded$seed, et.tsne.loaded@seed)
  expect_equal(et.tsne.loaded$selected_rows, et.tsne.loaded@selected_rows)
  expect_equal(et.tsne.loaded$selected_columns, et.tsne.loaded@selected_columns)
  expect_equal(et.tsne.loaded$tsne_result, et.tsne.loaded@tsne_result)
})

test_that("set input", {
  expect_equal(et.input.loaded$raw_counts, et.input.loaded@raw_counts)
  expect_equal(et.input.loaded$norm_counts, et.input.loaded@norm_counts)
  expect_equal(et.input.loaded$meta_data, et.input.loaded@meta_data)
  expect_equal(et.input.loaded$perplexity, et.input.loaded@perplexity)
  expect_equal(et.input.loaded$seed, et.input.loaded@seed)
  expect_equal(et.input.loaded$selected_rows, et.input.loaded@selected_rows)
  expect_equal(et.input.loaded$selected_columns, et.input.loaded@selected_columns)
})