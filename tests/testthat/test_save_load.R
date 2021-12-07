testthat::context("save and load")
library(expTSNE)
library(testthat)

test_counts = matrix(runif(10000), nrow = 10, ncol = 1000)
rownames(test_counts) = paste("row", seq_len(nrow(test_counts)))
colnames(test_counts) = paste("col", seq_len(ncol(test_counts)))

et.input = expTSNE.input(test_counts)
et.tsne = expTSNE.runTSNE(et.input)

expTSNE.save(et.input, "test_expTSNE.input", overwrite = TRUE)
expTSNE.save(et.tsne, "test_expTSNE", overwrite = TRUE)

et.input.loaded = expTSNE.load("test_expTSNE.input")
et.tsne.loaded = expTSNE.load("test_expTSNE")

test_that("save load equal", {
  expect_equal(colnames(et.input@raw_counts), colnames(et.input.loaded@raw_counts))
  expect_equal(rownames(et.input@raw_counts), rownames(et.input.loaded@raw_counts))
  expect_equal(et.input@raw_counts, et.input.loaded@raw_counts)
  
  expect_equal(colnames(et.input@norm_counts), colnames(et.input.loaded@norm_counts))
  expect_equal(rownames(et.input@norm_counts), rownames(et.input.loaded@norm_counts))
  expect_equal(et.input@norm_counts, et.input.loaded@norm_counts)
  
  expect_equal(colnames(et.input@meta_data), colnames(et.input.loaded@meta_data))
  expect_equal(rownames(et.input@meta_data), rownames(et.input.loaded@meta_data))
  expect_equal(et.input@meta_data, et.input.loaded@meta_data)
})

test_that("save load equal", {
  expect_equal(colnames(et.tsne@raw_counts), colnames(et.tsne.loaded@raw_counts))
  expect_equal(rownames(et.tsne@raw_counts), rownames(et.tsne.loaded@raw_counts))
  expect_equal(et.tsne@raw_counts, et.tsne.loaded@raw_counts)
  
  expect_equal(colnames(et.tsne@norm_counts), colnames(et.tsne.loaded@norm_counts))
  expect_equal(rownames(et.tsne@norm_counts), rownames(et.tsne.loaded@norm_counts))
  expect_equal(et.tsne@norm_counts, et.tsne.loaded@norm_counts)
  
  expect_equal(colnames(et.tsne@meta_data), colnames(et.tsne.loaded@meta_data))
  expect_equal(rownames(et.tsne@meta_data), rownames(et.tsne.loaded@meta_data))
  expect_equal(et.tsne@meta_data, et.tsne.loaded@meta_data)
})

unlink("test_expTSNE.input", recursive = TRUE)
unlink("test_expTSNE", recursive = TRUE)
