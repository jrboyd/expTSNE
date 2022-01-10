.names_expTSNE.input = function(x){
  c(
    "raw_counts",
    "norm_counts", 
    "meta_data", 
    "perplexity", 
    "seed", 
    "selected_rows", 
    "selected_columns"
  )
}
.names_expTSNE = function(x){
  c(.names_expTSNE.input(), "tsne_result")
}

### $ Accessor
setMethod("names", 
          "expTSNE.input",
          .names_expTSNE.input)
setMethod("names", 
          "expTSNE",
          .names_expTSNE)

.get_expTSNE = function(x, name){
  if(!name %in% names(x)) stop(name, " is not a valid item in this ", class(x), ".")
  switch (name,
          raw_counts = x@raw_counts,
          norm_counts = x@norm_counts, 
          meta_data = x@meta_data, 
          perplexity = x@perplexity, 
          seed = x@seed, 
          selected_rows = x@selected_rows, 
          selected_columns = x@selected_columns,
          tsne_result = x@tsne_result
  )
}

setMethod("$", 
          "expTSNE.input",
          .get_expTSNE
)

.set_expTSNE = function(x, name, value){
  if(!name %in% names(x)) stop(name, " is not a valid item in this ", class(x), ".")
  warn_msg = "This assignment is not supported.  No effect."
  switch (name,
          raw_counts = {
            x@raw_counts = value
          },
          norm_counts = {
            x@norm_counts = value
          }, 
          meta_data = {
            x@meta_data = value
          }, 
          perplexity = {
            x@perplexity = value
          }, 
          seed = {
            x@seed = value
          }, 
          selected_rows = {
            if(!all(value %in% rownames(x@norm_counts))) stop("Not all selected_rows were in rownames.")
            x@selected_rows = sort(unique(value))
          }, 
          selected_columns = {
            x@selected_columns = value
          },
          tsne_result = {
            x@tsne_result = value
          },
          {warning(warn_msg)}
          
  )
  x
}
setReplaceMethod("$", 
                 "expTSNE.input",
                 .set_expTSNE)
