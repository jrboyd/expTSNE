library(data.table)
library(expTSNE)



tab_files.target = dir("~/../dbgap/data/alignment_RNA-Seq/", pattern = ".ReadsPerGene.out.tab$", full.names = TRUE)
bam_files.target = dir("~/../dbgap/data/alignment_RNA-Seq/", pattern = "bam$", full.names = TRUE)

gtf_file = "~/indexes/HG38canon/GTF/gencode.v36.annotation.gtf"
sum(duplicated(substr(basename(bam_files), 1, 15)))


# strand_assessment = ssvRecipes::rnaseq_asses_strandedness(bam_files, gtf_file, max_bams = 20, sample_names = sub(".Aligned.sortedByCoord.out.bam", "", basename(bam_files)))
# strand_assessment.target = ssvRecipes::rnaseq_asses_strandedness("~/../dbgap/data/alignment_RNA-Seq/", gtf_file, max_bams = 20, sample_names = sub(".Aligned.sortedByCoord.out.bam", "", basename(bam_files)))
# 
# strand_assessment$scatter
# strand_assessment$tracks_plus
# strand_assessment$tracks_negative

sapply(tab_files.target[1:15], ssvRecipes::guess_lib_from_file)

mat.target = ssvRecipes::load_matrix_from_ReadsPerGene.out.tab(tab_files.target, lib_type = "unstranded")

dim(mat.target)

mat = mat.target
dim(mat)

# ssvQC::get_mapped_reads(c(bam_files, bam_files.target))

gene_max = apply(mat, 1, quantile, probs = .95)
hist(log10(gene_max+1))
k_expressed = gene_max > 100
sum(gene_max > 100)

mat.rpm = mat
for(i in seq_len(ncol(mat))){
  mat.rpm[, i] = mat[, i] /  sum(mat[,i]) * 1e6
}
# mat.rpm = apply(mat, 2, function(x){x/sum(x)*1e6})

dim(mat)
dim(mat.rpm)

mat = mat[k_expressed,]
mat.rpm = mat.rpm[k_expressed,]

dt = as.data.table(melt(mat))
setnames(dt, c("gene_id", "sample", "count"))
dt[, .(round(sum(count)/1e6, 1)), .(sample)][order(V1)]

et = expTSNE.input(raw_counts = mat, norm_counts = mat.rpm)

et

# et@norm_counts = (apply(et@raw_counts, 2, function(x)x/sum(x)*1e6))
dim(et@norm_counts)
dim(et@raw_counts)

parse_patient_ids = function(ids){
  dt = data.table(id = ids)
  dt = dt[, tstrsplit(id, "[\\.]", keep = 3)]
  dt = dt[, tstrsplit(V1, "-", keep = c(1:3))]
  dt[, paste(V1, V2, V3, sep = "-")]
}

parse_sample_codes = function(ids){
  dt = data.table(id = ids)
  dt = dt[, tstrsplit(id, "[\\.]", keep = 3)]
  dt = dt[, tstrsplit(V1, "-", keep = c(4))]
  sub("[A-Z]", "", dt$V1)
}

meta_dt = as.data.table(et$meta_data)
meta_dt[, patient_id := parse_patient_ids(column_id)]
meta_dt[, sample_code := parse_sample_codes(column_id)]

# 40,Recurrent Blood Derived Cancer - Peripheral Blood,TRB
# 03,Primary Blood Derived Cancer - Peripheral Blood,TB
# 04,Recurrent Blood Derived Cancer - Bone Marrow,TRBM
# 09,Primary Blood Derived Cancer - Bone Marrow,TBM
code2sample = c(
  "40" = "Recurrent Blood Derived Cancer - Peripheral Blood",
  "03" = "Primary Blood Derived Cancer - Peripheral Blood",
  "04" = "Recurrent Blood Derived Cancer - Bone Marrow",
  "09" = "Primary Blood Derived Cancer - Bone Marrow",
  "not_set" = "not_set"
)
meta_dt$sample_type = code2sample[meta_dt$sample_code]

target_dt = fread("~/R/SF_target_RNAseq/clinical_merged.all.csv")
clin_dt = target_dt[, .(patient_id, Gender, Race, Ethnicity, CNS.Status.at.Diagnosis, Cell.of.Origin)]

meta_dt =merge(meta_dt, clin_dt, by = "patient_id", all.x = TRUE)

msig.H = ssvRecipes::make_msigdb_TERM2GENE(species = "Homo sapiens", category = "H")
table(msig.H$gs_name)
sel_genes.gene_name = subset(msig.H, gs_name == "HALLMARK_COMPLEMENT")$gene_symbol

ref_gr = rtracklayer::import.gff(gtf_file, feature.type = "gene")
sel_genes = subset(ref_gr, gene_name %in% sel_genes.gene_name)$gene_id

sel_genes = intersect(sel_genes, rownames(et$norm_counts))
et$selected_rows = sel_genes
et$selected_rows
dim(et$norm_counts)
et = expTSNE.runTSNE(et)

et$meta_data = meta_dt

plot_expTSNE = function(et, color_var, facet_var = NULL, color_palette = "Dark2", show_all_points_per_facet = TRUE){
  pdt = merge(et$tsne_result, et$meta_data, by = "column_id")
  p = ggplot(pdt, aes_string(x = "tx", y = "ty", color = color_var)) 
  
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
  p = p +
    cowplot::theme_cowplot()
  p
}

ik_id = subset(ref_gr, gene_name == "IKZF1")$gene_id

et$meta_data$value = log10(et$norm_counts[ik_id,] + 1)
v = et$meta_data$value
v = (v-mean(v))/sd(v) 
v[v < -3] = -3
et$meta_data$value = v
et$meta_data


plot_expTSNE(et, color_var = "value")

plot_expTSNE(et, color_var = "Cell.of.Origin", facet_var = "Race")
plot_expTSNE(et, color_var = "value", facet_var = "Race") +
  labs(color = "IKZF1 expression")

plot_expTSNE(et, color_var = "value", facet_var = "Cell.of.Origin")

expTSNE.runApp(et)
