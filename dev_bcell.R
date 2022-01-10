library(data.table)
library(expTSNE)

counts <- read.csv("~/R/PR_Bcell_development/data/preB_healthy_counts.csv", header = TRUE, row.names = 1)

colnames(counts)
condition <-  c("proB", "preBI", "preBIIS", "preBIIL", "proB", "preBI", "preBIIS", "preBIIL", "proB", "preBI", "preBIIS", "preBIIL")
coldata <- data.frame(row.names=colnames(counts), condition)

tab_files = dir("~/R/PR_Bcell_development/data/", pattern = "tab$", full.names = TRUE)
bam_files = dir("~/R/PR_Bcell_development/data/", pattern = "bam$", full.names = TRUE)
tab_files.target = dir("~/../dbgap/data/alignment_RNA-Seq/", pattern = ".ReadsPerGene.out.tab$", full.names = TRUE)
bam_files.target = dir("~/../dbgap/data/alignment_RNA-Seq/", pattern = "bam$", full.names = TRUE)
tab_files.ccle = dir("~/R/SF_Ikaros_splicing/monaco_and_CCLE", pattern = ".ReadsPerGene.out.tab$", full.names = TRUE)
bam_files.ccle = dir("~/R/SF_Ikaros_splicing/monaco_and_CCLE", pattern = "bam$", full.names = TRUE)


gtf_file = "~/indexes/HG38canon/GTF/gencode.v36.annotation.gtf"
sum(duplicated(substr(basename(bam_files), 1, 15)))


# strand_assessment = ssvRecipes::rnaseq_asses_strandedness(bam_files, gtf_file, max_bams = 20, sample_names = sub(".Aligned.sortedByCoord.out.bam", "", basename(bam_files)))
# strand_assessment.target = ssvRecipes::rnaseq_asses_strandedness("~/../dbgap/data/alignment_RNA-Seq/", gtf_file, max_bams = 20, sample_names = sub(".Aligned.sortedByCoord.out.bam", "", basename(bam_files)))
# 
# strand_assessment$scatter
# strand_assessment$tracks_plus
# strand_assessment$tracks_negative

sapply(tab_files, ssvRecipes::guess_lib_from_file)
sapply(tab_files.target[1:15], ssvRecipes::guess_lib_from_file)

mat.dev = ssvRecipes::load_matrix_from_ReadsPerGene.out.tab(tab_files, lib_type = "unstranded")
write.table(mat.dev, file = "~/R/PR_Bcell_development/data/bcell_dev_and_careo.unstranded_counts.csv", sep = ",", quote = FALSE)
mat.target = ssvRecipes::load_matrix_from_ReadsPerGene.out.tab(tab_files.target, lib_type = "unstranded")

dim(mat.dev)
dim(mat.target)

common = intersect(rownames(mat.dev), rownames(mat.target))
length(common)

mat = cbind(mat.dev[common,], mat.target[common,])

dim(mat)

ssvQC::get_mapped_reads(bam_files, bam_files.target)

gene_max = apply(mat, 1, quantile, probs = .95)
hist(log10(gene_max+1))
k_expressed = gene_max > 100
sum(gene_max > 100)

mat.rpm = apply(mat, 2, function(x){x/sum(x)*1e6})

dim(mat)
dim(mat.rpm)

mat = mat[k_expressed,]
mat.rpm = mat.rpm[k_expressed,]

dt = as.data.table(melt(mat))
setnames(dt, c("gene_id", "sample", "count"))
dt[, .(round(sum(count)/1e6, 1)), .(sample)][order(V1)]

et = expTSNE.input(raw_counts = mat, norm_counts = mat.rpm)

et$meta_data$source = ifelse(grepl("donor", et$meta_data$column_id), "PR", "CAREO")
et$meta_data[grepl("TARGET", et$meta_data$column_id),]$source = "TARGET"

et

# et@norm_counts = (apply(et@raw_counts, 2, function(x)x/sum(x)*1e6))
dim(et@norm_counts)
dim(et@raw_counts)



table(et$meta_data$source)

subset(et$meta_data, source == "PR")
subset(et$meta_data, source == "CAREO")

msig.H = ssvRecipes::make_msigdb_TERM2GENE(species = "Homo sapiens", category = "H")
table(msig.H$gs_name)
sel_genes = subset(msig.H, gs_name == "HALLMARK_COMPLEMENT")$gene_symbol

et$selected_rows = sel_genes
et$selected_rows
dim(et$norm_counts)
et = expTSNE.runTSNE(et)

pdt = merge(et@tsne_result, et@meta_data, by = "column_id")
ggplot(pdt, aes(x = tx, y = ty, color = source)) + geom_point()


