testthat::context("run TSNE")
library(expTSNE)
library(testthat)

test_counts = matrix(runif(10000), nrow = 10, ncol = 1000)[, 1:150]
rownames(test_counts) = paste("row", seq_len(nrow(test_counts)))
colnames(test_counts) = paste("col", seq_len(ncol(test_counts)))

run_TSNE = expTSNE:::run_TSNE

test_that("run_tsne minimal", {
  tsne_df = run_TSNE(test_counts)
  expect_is(tsne_df, "data.frame")
})

test_that("run_tsne auto reduce perplexity", {
  expect_warning({run_TSNE(test_counts[,1:30])}, regexp = "auto reducing perplexity")
})

test_that("run_tsne normalization", {
  test_counts[, 1:5] = test_counts[, 1:5]*10
  tsne_df.norm = run_TSNE(test_counts, apply_normalization = TRUE)
  tsne_df.norm$spike = c(rep("spike", 5), rep("nope", 145))
  ggplot(tsne_df.norm, aes(x = tx, y = ty, color = spike)) + geom_point()
  expect_is(tsne_df.norm, "data.frame")
  
  
  if(FALSE){
    #without normalization the spikes all cluster together
    tsne_df = run_TSNE(test_counts, apply_normalization = FALSE)
    tsne_df$spike = c(rep("spike", 5), rep("nope", 145))
    ggplot(tsne_df, aes(x = tx, y = ty, color = spike)) + geom_point()
  }
})
