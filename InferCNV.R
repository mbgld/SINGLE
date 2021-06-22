# Set-up
library(Seurat); library(infercnv)
Infer_CNV.outdir = "Results/InferCNV/"

### Part	A. Set input files
### 1 Create read count matrix file from Seurat.obj
Infercnv.matrix.file <- paste0(Infer_CNV.outdir, "inferCNV.matrix")
write.table(Case_{}.sfj@assays$RNA@counts, file = Infercnv.matrix.file, quote = F, sep = "\t")

### 2 Create ann.file from Seurat.obj
Infercnv.ann.file <- paste0(Infer_CNV.outdir, "inferCNV.annotation")
infercnv.ann.df <- data.frame(barcode = rownames(Case_{}.sfj@meta.data), clusters = Case_{}.sfj@active.ident)
write.table(x = infercnv.ann.df, file = Infercnv.ann.file, quote = F, sep = "\t", row.names = F, col.names = F) 

### Part B. Run inferCNV
# 1.Indicate gtf file location
gene.order.file <- "/media/yschang/T/other_materials/hg19.inferCNV.gtf"

# 2 Check DimPlot and indicate ref.group
DimPlot(Case_{}.sfj, reduction = "umap", label = T)
levels(Case_{}.sfj@active.ident)
ref.group <- c("NL_Alveolar_cells", "0.T-cells", "1.T-cells", "3.Dendritic cells", "4.Fibroblasts", "5.Endothelial cells", "6.Fibroblasts", "7.Bronchial epithelial cells", "8.Mast cells") # It should be adjusted from individual case.

# 3. Create infer cnv objects
infercnv.obj <- CreateInfercnvObject(raw_counts_matrix = paste0(Infer_CNV.outdir,'inferCNV.matrix'), gene_order_file = gene.order.file, annotations_file = paste0(Infer_CNV.outdir,'inferCNV.annotation'), ref_group_names = ref.group, delim = "\t", chr_exclude = c("X", "Y"))

# 4. run infercnv with fine tuned parameters; analysis_mode
infercnv_run = infercnv::run(infercnv_obj = infercnv.obj, num_threads = 24, cutoff=0.1, out_dir= Infer_CNV.outdir, analysis_mode="subclusters", cluster_by_groups=F, plot_steps=T, denoise=T, no_prelim_plot=F, k_obs_groups = 3, HMM=T)

# 5. Plot infercnv
plot_cnv(infercnv_obj = infercnv.obj, out_dir = New_infer_CNV.outdir, title = "InferCNV_Plot", obs_title = "Lung cancer", ref_title = "Normal cells", output_format = "png", png_res = 1200)

### Part C. save
save.image(paste0(Infer_CNV.outdir, "Case_{}_InferCNV.RData"))
gc();q()
