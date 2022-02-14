#' Title
#'
#' @param et expTSNE object
#' @param color_var 
#' @param facet_var 
#' @param color_palette 
#' @param show_all_points_per_facet 
#'
#' @return
#' @export
#'
#' @examples
#' ex_data = system.file("extdata/test_expTSNE", package = "expTSNE", mustWork = TRUE)
#' et = expTSNE.load(ex_data)
#' plot_expTSNE(et)
plot_expTSNE = function(et, color_var = NULL, facet_var = NULL, color_palette = "Dark2", show_all_points_per_facet = TRUE){
  et = prep_expTSNE_for_gene(et, color_var)
  et = prep_expTSNE_for_gene(et, facet_var, bins = 4)
  pdt = merge(et$tsne_result, et$meta_data, by = "column_id")
  p = ggplot(pdt, aes_string(x = "tx", y = "ty", color = color_var, id = "column_id")) 
  
  if(is.null(color_var)){
    if(!is.null(facet_var)){
      if(show_all_points_per_facet){
        p = p +
          annotate("point", x = pdt$tx, y = pdt$ty, color = "gray", size = .5)
      }
      p = p + 
        geom_point() + 
        facet_wrap(paste0("~", facet_var))
      
    }else{
      p = p + 
        geom_point() 
    }
  }else{
    if(is.numeric(pdt[[color_var]])){
      sc = scale_color_viridis_c() 
    }else{
      sc = scale_color_brewer(palette = color_palette) 
    }
    if(!is.null(facet_var)){
      if(show_all_points_per_facet){
        p = p +
          annotate("point", x = pdt$tx, y = pdt$ty, color = "gray", size = .5)
      }
      p = p + 
        geom_point() + 
        sc + 
        facet_wrap(paste0("~", facet_var))
      
    }else{
      p = p + 
        geom_point() + 
        sc
    }
  }
  p = p +
    cowplot::theme_cowplot()
  p
}

prep_expTSNE_for_gene = function(et, var, bins = 0){
  if(is.null(var)) return(et)
  if(!var %in% colnames(et$meta_data) & var %in% rownames(et$norm_counts)){
    vals = et$norm_counts[var,]
    if(bins < 1){
      et$meta_data[[var]] = vals[et$meta_data$column_id]   
    }else{
      brks = quantile(vals, seq(0, bins)/bins)
      et$meta_data[[var]] = cut(vals[et$meta_data$column_id], breaks = brks, include.lowest = TRUE)
    }
  }
  et
}
