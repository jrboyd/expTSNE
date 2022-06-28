library(TARGETdata)
library(expTSNE)
options(mc.cores = 10)
TARGETdata.runApp.preload()


mrna_rpm_mat.full[1:5, 1:5]

meta_dt = make_meta_dt(mrna_rpm_mat.full,  clin_dt)
meta_dt = as.data.frame(meta_dt)
rownames(meta_dt) = meta_dt$sample_id

mrna_rpm_mat.full = filter_expression_to_valid(mrna_rpm_mat.full, meta_dt)
mrna_count_mat.full = filter_expression_to_valid(mrna_count_mat.full, meta_dt)



et = expTSNE.input(mrna_count_mat.full, log2(mrna_rpm_mat.full + .1), meta_data = meta_dt, norm_description = "log2 of RPM")
et = expTSNE.runTSNE(et)

expTSNE.runApp(et)

# et@meta_data$column_id = et@meta_data$sample_id

expTSNE::plot_expTSNE(et, color_var = "cluster_id", facet_var = "ik_status")
