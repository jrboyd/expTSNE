#one time modify to example data to provide grouping to data
if(FALSE){
  ex_data = system.file("extdata/test_expTSNE.input", package = "expTSNE", mustWork = TRUE)
  et = expTSNE.load(ex_data)
  ids = colnames(et$norm_counts)
  set.seed(0)
  grp_a = sample(ids, 300)
  grp_b = sample(ids, 400)
  grp_b = setdiff(grp_b, grp_a)
  
  et$meta_data$group = "C"
  et$meta_data[rownames(et$meta_data) %in% grp_b,]$group = "B"
  et$meta_data[rownames(et$meta_data) %in% grp_a,]$group = "A"
  table(et$meta_data$group )
  
  et$raw_counts[, colnames(et$raw_counts) %in% grp_b][1:3,] = et$raw_counts[, colnames(et$raw_counts) %in% grp_b][1:3,] *1.5
  et$raw_counts[, colnames(et$raw_counts) %in% grp_a][3:4,] = et$raw_counts[, colnames(et$raw_counts) %in% grp_a][3:4,] *1.5
  et$norm_counts = et$raw_counts
  expTSNE.save(et, "inst/extdata/test_expTSNE.input/", overwrite = TRUE)
  
  et = expTSNE.runTSNE(et)
  expTSNE.save(et, "inst/extdata/test_expTSNE", overwrite = TRUE)
}