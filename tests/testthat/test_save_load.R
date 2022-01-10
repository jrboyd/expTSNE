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

test_that("save load equal for expTSNE.input", {
  expect_equal(class(et.input), class(et.input.loaded))
  
  expect_equal(colnames(et.input@raw_counts), colnames(et.input.loaded@raw_counts))
  expect_equal(rownames(et.input@raw_counts), rownames(et.input.loaded@raw_counts))
  expect_equal(et.input@raw_counts, et.input.loaded@raw_counts)
  
  expect_equal(colnames(et.input@norm_counts), colnames(et.input.loaded@norm_counts))
  expect_equal(rownames(et.input@norm_counts), rownames(et.input.loaded@norm_counts))
  expect_equal(et.input@norm_counts, et.input.loaded@norm_counts)
  
  expect_equal(colnames(et.input@meta_data), colnames(et.input.loaded@meta_data))
  expect_equal(rownames(et.input@meta_data), rownames(et.input.loaded@meta_data))
  expect_equal(et.input@meta_data, et.input.loaded@meta_data)
  
  expect_equal(et.input@perplexity, et.input.loaded@perplexity)
  
  expect_equal(et.input@seed, et.input.loaded@seed)
  
  expect_equal(et.input@selected_rows, et.input.loaded@selected_rows)
  
  expect_equal(et.input@selected_columns, et.input.loaded@selected_columns)
})

test_that("save load equal for expTSNE", {
  expect_equal(class(et.tsne), class(et.tsne.loaded))
  
  expect_equal(colnames(et.tsne@raw_counts), colnames(et.tsne.loaded@raw_counts))
  expect_equal(rownames(et.tsne@raw_counts), rownames(et.tsne.loaded@raw_counts))
  expect_equal(et.tsne@raw_counts, et.tsne.loaded@raw_counts)
  
  expect_equal(colnames(et.tsne@norm_counts), colnames(et.tsne.loaded@norm_counts))
  expect_equal(rownames(et.tsne@norm_counts), rownames(et.tsne.loaded@norm_counts))
  expect_equal(et.tsne@norm_counts, et.tsne.loaded@norm_counts)
  
  expect_equal(colnames(et.tsne@meta_data), colnames(et.tsne.loaded@meta_data))
  expect_equal(rownames(et.tsne@meta_data), rownames(et.tsne.loaded@meta_data))
  expect_equal(et.tsne@meta_data, et.tsne.loaded@meta_data)
  
  expect_equal(et.tsne@perplexity, et.tsne.loaded@perplexity)
  
  expect_equal(et.tsne@seed, et.tsne.loaded@seed)
  
  expect_equal(et.tsne@selected_rows, et.tsne.loaded@selected_rows)
  
  expect_equal(et.tsne@selected_columns, et.tsne.loaded@selected_columns)
  
  expect_equal(colnames(et.tsne@tsne_result), colnames(et.tsne.loaded@tsne_result))
  expect_equal(rownames(et.tsne@tsne_result), rownames(et.tsne.loaded@tsne_result))
  expect_equal(et.tsne@tsne_result$tx, et.tsne.loaded@tsne_result$tx)
  expect_equal(et.tsne@tsne_result$ty, et.tsne.loaded@tsne_result$ty)
  expect_equal(et.tsne@tsne_result$column_id, et.tsne.loaded@tsne_result$column_id)
})

unlink("test_expTSNE.input", recursive = TRUE)
unlink("test_expTSNE", recursive = TRUE)
