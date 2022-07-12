library(TARGETdata)
library(expTSNE)
options(mc.cores = 10)
TARGETdata.runApp.preload()

if(!exists("et")){
  mrna_rpm_mat.full[1:5, 1:5]
  
  meta_dt = make_meta_dt(mrna_rpm_mat.full,  clin_dt)
  meta_dt = as.data.frame(meta_dt)
  rownames(meta_dt) = meta_dt$sample_id
  
  mrna_rpm_mat.full = filter_expression_to_valid(mrna_rpm_mat.full, meta_dt)
  mrna_count_mat.full = filter_expression_to_valid(mrna_count_mat.full, meta_dt)
  
  
  
  et = expTSNE.input(mrna_count_mat.full, log2(mrna_rpm_mat.full + .1), meta_data = meta_dt, norm_description = "log2 of RPM")
  et = expTSNE.runTSNE(et)
} 
expTSNE.runApp(et, display.mode="showcase")

# et@meta_data$column_id = et@meta_data$sample_id

expTSNE::plot_expTSNE(et, color_var = "cluster_id", facet_var = "ik_status")
nn_clust = expTSNE:::nn_clust

runApp(list(
  ui = fluidPage(
    theme = "bootstrap.css",
    # Application title
    shinyjs::useShinyjs(),
    titlePanel("TCGA t-sne"),
    expTSNE:::ui_upload(id = "up1"),
    expTSNE:::ui_point_selection(n_items = ncol(et$raw_counts)),
    expTSNE:::ui_goi()
    # expTSNE:::ui_upload(id = "up2")
  ), 
  
  
  server = function(input, output, session) {
    
    #all added custom gene sets
    custom_gene_sets = reactiveVal(list())
    expression_data = reactiveVal(et$norm_counts)
    
    expTSNE:::server_upload(
      id = "up1",
      expression_data = expression_data
    )
    
    # tsne_clust, meta_data, tsne_res
    meta_data = reactiveVal(et$meta_data)
    tsne_res = reactiveVal(as.data.table(et$tsne_result))
    tsne_clust = reactiveVal()
    # expTSNE:::server_point_selection(tsne_res = tsne_res, meta_data = meta_data, tsne_clust = tsne_clust)
    expTSNE:::server_point_selection(tsne_res = tsne_res, meta_data = meta_data, tsne_clust = tsne_clust)
    expTSNE:::server_goi()
    # expTSNE:::server_upload(
    #   id = "up2",
    #   expression_data = expression_data
    # )
  }
))
