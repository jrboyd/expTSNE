setClassUnion("numericOrNull", c("numeric", "NULL"))

.check_valid_expTSNE = function(object){
  msg = character()
  if(is.null(rownames(object@raw_counts))) msg = c(msg, "raw_counts must have rownames.")
  if(is.null(colnames(object@raw_counts))) msg = c(msg, "raw_counts must have colnames.")
  if(is.null(object@selected_rows)) msg = c(msg, "selected_rows cannot be NULL.")
  if(is.null(object@selected_columns)) msg = c(msg, "selected_columns cannot be NULL.")
  if(!all(object@selected_rows %in% rownames(object@raw_counts))) msg = c(msg, "selected_rows must all be in rownames.")
  if(!all(object@selected_columns %in% colnames(object@raw_counts))) msg = c(msg, "selected_columns must all be in colnames.")
  if(!all(ncol(object@raw_counts) == ncol(object@norm_counts))) msg = c(msg, "norm_counts and raw_counts colnames must be equal.")
  if(!all(nrow(object@raw_counts) == nrow(object@norm_counts))) msg = c(msg, "norm_counts and raw_counts rownames must be equal.")
  if(!all(colnames(object@raw_counts) == colnames(object@norm_counts))) msg = c(msg, "norm_counts and raw_counts colnames must be equal.")
  if(!all(rownames(object@raw_counts) == rownames(object@norm_counts))) msg = c(msg, "norm_counts and raw_counts rownames must be equal.")
  if(nrow(object@meta_data) != ncol(object@raw_counts)) msg = c(msg, "meta_data must have one row per column in raw_counts.")
  if(is.null(rownames(object@meta_data))) msg = c(msg, "meta_data must have rownames.")
  if(!all(rownames(object@meta_data) == colnames(object@raw_counts))) msg = c(msg, "meta_data rownames must equal raw_counts colnames.")
  if(length(msg) > 0){
    return(msg)
  }else{
    return(TRUE)
  }
}

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
           norm_description = "character",
           meta_data = "data.frame",
           column_id_var = "character",
           perplexity = "numeric",
           seed = "numericOrNull",
           selected_rows = "character",
           selected_columns = "character"
         ),
         validity = .check_valid_expTSNE
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
    norm_description = "normalized counts",
    meta_data = data.frame(column_id = colnames(raw_counts), row.names = colnames(raw_counts)),
    column_id_var = "column_id",
    perplexity = 30,
    seed = 0,
    selected_rows = rownames(raw_counts),
    selected_columns = colnames(raw_counts)
){
  if(is.null(rownames(raw_counts))) stop("raw_counts must have rownames.")
  if(is.null(colnames(raw_counts))) stop("raw_counts must have colnames.")
  if(!all(ncol(raw_counts) == ncol(norm_counts))) stop("norm_counts and raw_counts colnames must be equal.")
  if(!all(nrow(raw_counts) == nrow(norm_counts))) stop("norm_counts and raw_counts rownames must be equal.")
  if(!all(colnames(raw_counts) == colnames(norm_counts))) stop("norm_counts and raw_counts colnames must be equal.")
  if(!all(rownames(raw_counts) == rownames(norm_counts))) stop("norm_counts and raw_counts rownames must be equal.")
  if(is.data.table(meta_data)) meta_data = as.data.frame(meta_data)
  if(is.null(rownames(meta_data))){
    if(is.null(meta_data[[column_id_var]])){
      stop(column_id_var, " must be in meta_data if rownames not set")
    }
    rownames(meta_data) = meta_data[[column_id_var]]
  }
  if(is.null(meta_data[[column_id_var]])){
    if(is.null(rownames(meta_data))){
      stop(column_id_var, " must be in meta_data if rownames not set")
    }
    meta_data[[column_id_var]] = rownames(meta_data)
  }
  new("expTSNE.input", 
      raw_counts = raw_counts,
      norm_counts = norm_counts,
      norm_description = norm_description,
      meta_data = meta_data,
      perplexity = perplexity,
      seed = seed,
      selected_rows = selected_rows,
      selected_columns = selected_columns,
      column_id_var = column_id_var)
}

#' expTSNE
#'
#' @slot tsne_result data.frame. 
#'
#' @return
#' @export
#' @rdname expTSNE
setClass("expTSNE", 
         slots = c(
           tsne_result = "data.frame"
         ),
         contains = "expTSNE.input"
)

#' get_args
#'
#' returns parameters of calling function as a named list.
#'
#' @param env
#' @param ...
#'
#' @return
#'
#' @examples
get_args = function(env = parent.frame(), to_ignore = character(), ...){
  args = c(as.list(env), list(...))
  args = args[!names(args) %in% to_ignore]
  args[order(names(args))]
}
#' digest_args
#'
#' returns digest results of name list of parameters of calling function
#'
#' @param env
#' @param ...
#'
#' @return
#' @import digest
#'
#' @examples
digest_args = function(env = parent.frame(), to_ignore = character(), ...){
  digest::digest(get_args(env, to_ignore, ...))
}


#' @param bfc
#' @param rname
#' @param FUN
#' @param version
#' @param force_overwrite
#'
#' @return
#' @rawNamespace import(BiocFileCache, except = show)
#'
#' @examples
bfcif = function(bfc, rname, FUN,
                 version = "v1",
                 force_overwrite = getOption("TSNE_FORCE_CACHE_OVERWRITE", FALSE)){
  # is rname in cache?
  vrname = paste0(rname, "_", version)
  if(nrow(BiocFileCache::bfcquery(bfc, query = vrname, field = "rname")) == 0){
    cache_path = BiocFileCache::bfcnew(bfc, rname = vrname)
    
  }else{
    cache_path = BiocFileCache::bfcrpath(bfc, vrname)
  }
  # does cached file exist?
  if(file.exists(cache_path) && !force_overwrite){
    load(BiocFileCache::bfcrpath(bfc, vrname))
  }else{
    res = FUN()
    save(res, file = cache_path)
  }
  # return either new results or cached results
  res
}


#' run_TSNE
#'
#' @param counts 
#' @param perplexity 
#' @param seed 
#'
#' @return
#' @import Rtsne
#' @export
#' @examples
#' ex_data = system.file("extdata/test_expTSNE.input", package = "expTSNE", mustWork = TRUE)
#' et = expTSNE.load(ex_data)
#' tsne_df = run_TSNE(et$norm_counts)
run_TSNE = function(counts, 
                    apply_normalization = FALSE, 
                    perplexity = 30, 
                    seed = NULL, 
                    column_id_var = "column_id",
                    bfc = getOption("TSNE_BFC", BiocFileCache::BiocFileCache())){
  rname = digest_args(to_ignore = c("bfc", "column_id"))
  
  .run_TSNE_FUN = function(){
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
    tsne_res
  }
  tsne_res = bfcif(bfc = bfc, 
                   rname = rname, 
                   FUN = .run_TSNE_FUN)
  
  
  
  tsne_df = as.data.table(tsne_res$Y)
  colnames(tsne_df) = c("tx", "ty")
  tsne_df[[column_id_var]] = colnames(counts)
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
#' ex_data = system.file("extdata/test_expTSNE.input", package = "expTSNE", mustWork = TRUE)
#' et = expTSNE.load(ex_data)
#' et.ran = expTSNE.runTSNE(et)
expTSNE.runTSNE = function(et){
  tsne_df = run_TSNE(et@norm_counts[et$selected_rows, et$selected_columns], perplexity = et$perplexity, seed = et$seed, column_id_var = et$column_id_var)
  new("expTSNE", 
      tsne_result = tsne_df,
      raw_counts = et$raw_counts,
      norm_counts = et$norm_counts,
      norm_description = et$norm_description,
      meta_data = et$meta_data,
      perplexity = et$perplexity,
      seed = et$seed,
      selected_rows = et$selected_rows,
      selected_columns = et$selected_columns,
      column_id_var = et$column_id_var)
}

#' expTSNE.save
#'
#' @param et expTSNE or expTSNE.input object to save
#' @param save_dir directory to save to. Will not be overwritten by default if it exists.
#' @param overwrite if TRUE, overwrite contents of save_dir.
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
      if(!all(dir(save_dir) %in% c(
        "meta_data.csv",
        "norm_counts.csv",
        "perplexity.txt",
        "raw_counts.csv",
        "seed.txt",
        "selected_columns.txt",
        "selected_rows.txt",
        "tsne_result.csv"
      ))){
        stop("Output contains non-standard files and cannot be overwritten safely. Please specify new save location.")
      }else{
        unlink(save_dir, recursive = TRUE)  
      }
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

.show_expTSNE.input = function(object){
  message("I am an expTSNE.input")
}
setMethod("show", "expTSNE.input", .show_expTSNE.input)

.show_expTSNE = function(object){
  message("I am an expTSNE")
}
setMethod("show", "expTSNE", .show_expTSNE)