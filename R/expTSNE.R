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
#' @rdname expTSNE
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
  tsne_df$column_id = colnames(counts)
  rownames(tsne_df) = colnames(counts)
  tsne_df
}

#' expTSNE.runTSNE
#'
#' @param et 
#'
#' @return
#' @export
#' @rdname expTSNE
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

# setGeneric("expTSNE.save")
# 
# setMethod("expTSNE.save", signature = c("expTSNE.input", "character", "logical"), function(et, save_dir, overwrite){
#   if(dir.exists(save_dir)){
#     if(!overwrite){
#       stop("Output already exists.  Run with `overwite = TRUE` to replace.")  
#     }else{
#       unlink(save_dir, recursive = TRUE)
#     }
#   }
#   dir.create(save_dir, recursive = TRUE)
#   write.table(et@raw_counts, file.path(save_dir, "raw_counts.csv"), sep = ",", quote = FALSE)
#   write.table(et@norm_counts, file.path(save_dir, "norm_counts.csv"), sep = ",", quote = FALSE)
#   write.table(et@meta_data, file.path(save_dir, "meta_data.csv"), sep = ",", quote = FALSE)
#   write(et@perplexity, file = file.path(save_dir, "perplexity.txt"))
#   if(is.null(et@seed)){
#     write("NULL", file = file.path(save_dir, "seed.txt"))  
#   }else{
#     write(et@seed, file = file.path(save_dir, "seed.txt"))  
#   }
#   write(et@selected_rows, file = file.path(save_dir, "selected_rows.txt"))
#   write(et@selected_columns, file = file.path(save_dir, "selected_columns.txt"))
# })
# 
# setMethod("expTSNE.save", signature = c("expTSNE", "character", "logical"), function(et, save_dir, overwrite){
#   if(dir.exists(save_dir)){
#     if(!overwrite){
#       stop("Output already exists.  Run with `overwite = TRUE` to replace.")  
#     }else{
#       unlink(save_dir, recursive = TRUE)
#     }
#   }
#   dir.create(save_dir, recursive = TRUE)
#   write.table(et@raw_counts, file.path(save_dir, "raw_counts.csv"), sep = ",", quote = FALSE)
#   write.table(et@norm_counts, file.path(save_dir, "norm_counts.csv"), sep = ",", quote = FALSE)
#   write.table(et@meta_data, file.path(save_dir, "meta_data.csv"), sep = ",", quote = FALSE)
#   write(et@perplexity, file = file.path(save_dir, "perplexity.txt"))
#   if(is.null(et@seed)){
#     write("NULL", file = file.path(save_dir, "seed.txt"))  
#   }else{
#     write(et@seed, file = file.path(save_dir, "seed.txt"))  
#   }
#   write(et@selected_rows, file = file.path(save_dir, "selected_rows.txt"))
#   write(et@selected_columns, file = file.path(save_dir, "selected_columns.txt"))
#   write(et@selected_columns, file = file.path(save_dir, "selected_columns.txt"))
# })

#' expTSNE.save
#'
#' @param et
#' @param save_dir
#' @param overwrite
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
#' et.run = expTSNE.runTSNE(et)
#'
#' expTSNE.save(et, "test_output")
#' expTSNE.save(et.run, "test_output.tsne")
expTSNE.save = function(et, save_dir, overwrite = FALSE){
  if(dir.exists(save_dir)){
    if(!overwrite){
      stop("Output already exists.  Run with `overwite = TRUE` to replace.")
    }else{
      unlink(save_dir, recursive = TRUE)
    }
  }
  dir.create(save_dir, recursive = TRUE)
  write.table(et@raw_counts, file.path(save_dir, "raw_counts.csv"), sep = ",", quote = FALSE)
  write.table(et@norm_counts, file.path(save_dir, "norm_counts.csv"), sep = ",", quote = FALSE)
  write.table(et@meta_data, file.path(save_dir, "meta_data.csv"), sep = ",", quote = FALSE)
  write(et@perplexity, file = file.path(save_dir, "perplexity.txt"))
  if(is.null(et@seed)){
    write("NULL", file = file.path(save_dir, "seed.txt"))
  }else{
    write(et@seed, file = file.path(save_dir, "seed.txt"))
  }
  write(et@selected_rows, file = file.path(save_dir, "selected_rows.txt"))
  write(et@selected_columns, file = file.path(save_dir, "selected_columns.txt"))
  
  if(is(et, "expTSNE")){
    write.table(et@tsne_result, file.path(save_dir, "tsne_result.csv"), sep = ",", quote = FALSE)  
  }
}

#' expTSNE.load
#'
#' @param et 
#' @param save_dir 
#'
#' @return
#' @export
#' @rdname expTSNE
#' @examples
#' et = expTSNE.load("test_output")
#' et.tsne = expTSNE.load("test_output.tsne")
expTSNE.load = function(save_dir){
  raw_counts = as.matrix(read.table(file.path(save_dir, "raw_counts.csv"), sep = ",", as.is = TRUE, check.names = FALSE))
  norm_counts = as.matrix(read.table(file.path(save_dir, "norm_counts.csv"), sep = ",", as.is = TRUE, check.names = FALSE))
  meta_data = read.table(file.path(save_dir, "meta_data.csv"), sep = ",", check.names = FALSE)
  perplexity = read.table(file.path(save_dir, "perplexity.txt"))[[1]]
  seed = read.table(file.path(save_dir, "seed.txt"), as.is = TRUE)[[1]]
  if(seed == "NULL") seed = NULL
  selected_rows = read.table(file.path(save_dir, "selected_rows.txt"), as.is = TRUE, sep = "\t")[[1]]
  selected_columns = read.table(file.path(save_dir, "selected_columns.txt"), as.is = TRUE, sep = "\t")[[1]]
  
  if(file.exists(file.path(save_dir, "tsne_result.csv"))){
    tsne_result = read.table(file.path(save_dir, "tsne_result.csv"), sep = ",", as.is = TRUE, check.names = FALSE)
    new("expTSNE",
      raw_counts = raw_counts,
      norm_counts = norm_counts, 
      meta_data = meta_data, 
      perplexity = perplexity, 
      seed = seed, 
      selected_rows = selected_rows, 
      selected_columns = selected_columns,
      tsne_result = tsne_result
    )
  }else{
    expTSNE.input(
      raw_counts = raw_counts,
      norm_counts = norm_counts, 
      meta_data = meta_data, 
      perplexity = perplexity, 
      seed = seed, 
      selected_rows = selected_rows, 
      selected_columns = selected_columns
    )
  }
}

setMethod("show", "expTSNE.input", function(object){
  message("I am an expTSNE.input")
})

setMethod("show", "expTSNE", function(object){
  message("I am an expTSNE")
})