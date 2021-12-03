#https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/Extensions.html
# https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/SummarizedExperiment/inst/doc/Extensions.html
#' expTSNE
#'
#' @slot raw_counts matrix. 
#' @slot norm_counts matrix. 
#' @slot meta_data data.frame. 
#' @slot perplexity numeric. 
#' @slot seed numeric. 
#' @slot selected_rows character. 
#' @slot selected_columns character. 
#'
#' @return
#' @export
#' @rdname expTSNE
setClass("expTSNE", 
         slots = c(
           raw_counts = "matrix",
           norm_counts = "matrix",
           meta_data = "data.frame",
           perplexity = "numeric",
           seed = "numeric",
           selected_rows = "character",
           selected_columns = "character"
           ))

#' expTSNE
#'
#' @return
#' @export
#' @rdname expTSNE
#' @examples
expTSNE = function(
  raw_counts,
  norm_counts = raw_counts,
  meta_data = data.frame(names = colnames(raw_counts)),
  perplexity = 30,
  seed = 0,
  selected_rows = rownames(raw_counts),
  selected_columns = colnames(raw_counts)
){
  new("expTSNE", 
      raw_counts = raw_counts,
      norm_counts = norm_counts,
      meta_data = meta_data,
      perplexity = perplexity,
      seed = seed,
      selected_rows = selected_rows,
      selected_columns = selected_columns)
}


