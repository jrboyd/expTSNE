setClassUnion("numericOrNull", c("numeric", "NULL"))

#https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/Extensions.html
# https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/SummarizedExperiment/inst/doc/Extensions.html
#' expTSNE.input
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
setClass("expTSNE.input", 
         slots = c(
           raw_counts = "matrix",
           norm_counts = "matrix",
           meta_data = "data.frame",
           perplexity = "numeric",
           seed = "numericOrNull",
           selected_rows = "character",
           selected_columns = "character"
         )
)

#' expTSNE.input
#'
#' @return
#' @export
#' @rdname expTSNE
#' @examples
#' test_counts = matrix(runif(10000), nrow = 10, ncol = 1000)
#' rownames(test_counts) = paste("row", seq_len(nrow(test_counts)))
#' colnames(test_counts) = paste("col", seq_len(ncol(test_counts)))
#' 
#' et = expTSNE.input(test_counts)
#' et = expTSNE.runTSNE(et)
#' expTSNE.runApp(et)
expTSNE.input = function(
  raw_counts,
  norm_counts = raw_counts,
  meta_data = data.frame(column_id = colnames(raw_counts), row.names = colnames(raw_counts)),
  perplexity = 30,
  seed = 0,
  selected_rows = rownames(raw_counts),
  selected_columns = colnames(raw_counts)
){
  if(is.null(rownames(raw_counts))) stop("raw_counts must have rownames.")
  if(is.null(colnames(raw_counts))) stop("raw_counts must have colnames.")
  if(is.null(selected_rows)) stop("selected_rows cannot be NULL.")
  if(is.null(selected_columns)) stop("selected_columns cannot be NULL.")
  if(!all(selected_rows %in% rownames(raw_counts))) stop("selected_rows must all be in rownames.")
  if(!all(selected_columns %in% colnames(raw_counts))) stop("selected_columns must all be in colnames.")
  if(!all(ncol(raw_counts) == ncol(norm_counts))) stop("norm_counts and raw_counts colnames must be equal.")
  if(!all(nrow(raw_counts) == nrow(norm_counts))) stop("norm_counts and raw_counts rownames must be equal.")
  if(!all(colnames(raw_counts) == colnames(norm_counts))) stop("norm_counts and raw_counts colnames must be equal.")
  if(!all(rownames(raw_counts) == rownames(norm_counts))) stop("norm_counts and raw_counts rownames must be equal.")
  if(nrow(meta_data) != ncol(raw_counts)) stop("meta_data must have one row per column in raw_counts.")
  if(is.null(rownames(meta_data))) stop("meta_data must have rownames.")
  if(!all(rownames(meta_data) == colnames(raw_counts))) stop("meta_data rownames must equal raw_counts colnames.")
  new("expTSNE.input", 
      raw_counts = raw_counts,
      norm_counts = norm_counts,
      meta_data = meta_data,
      perplexity = perplexity,
      seed = seed,
      selected_rows = selected_rows,
      selected_columns = selected_columns)
}

#' expTSNE
#'
#' @slot tsne_result data.frame. 
#'
#' @return
#' @export
#'
#' @examples
setClass("expTSNE", 
         slots = c(
           tsne_result = "data.frame"
         ),
         contains = "expTSNE.input"
)

#' run_TSNE
#'
#' @param counts 
#' @param perplexity 
#' @param seed 
#'
#' @return
#' @import Rtsne
#' @examples
run_TSNE = function(counts, apply_normalization = FALSE, perplexity = 30, seed = NULL){
  set.seed(seed)
  if(apply_normalization){
    counts = Rtsne::normalize_input(counts)
  }
  if(perplexity > ncol(counts)/4){
    warning("auto reducing perplexity")
    perplexity = round(ncol(counts)/4)
  }
  tsne_res = Rtsne::Rtsne(t(counts), perplexity = perplexity, num_threads = getOption("mc.cores", 1), check_duplicates = FALSE)
  set.seed(NULL)

  tsne_df = as.data.table(tsne_res$Y)
  colnames(tsne_df) = c("tx", "ty")
  tsne_df$sample_id = colnames(counts)
  tsne_df
}

#' expTSNE.runTSNE
#'
#' @param et 
#'
#' @return
#' @export
#'
#' @examples
expTSNE.runTSNE = function(et){
  tsne_df = run_TSNE(et@norm_counts, perplexity = et@perplexity, seed = et@seed)
  new("expTSNE", 
      tsne_result = tsne_df,
      raw_counts = et@raw_counts,
      norm_counts = et@norm_counts,
      meta_data = et@meta_data,
      perplexity = et@perplexity,
      seed = et@seed,
      selected_rows = et@selected_rows,
      selected_columns = et@selected_columns)
}