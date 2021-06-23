# set-up
library(Seurat); library(dplyr); library(multtest);library(metap); library(monocle); library(plyr)
Trj.outdir = './Results/Trajectory/'

### Part A. Assign id
# 1. Keep subcluster.id
TC.sfj$subcluster.id <- Idents(TC.sfj)
TC.sfj <- RenameIdents(TC.sfj, '0' = 'naive CD4', '1' = 'naive CD8', '2' = 'naive CD8', '3'= 'NKT', '4' = 'CD4_memory', '5' = 'CD8_effector', '6' = 'CD8_effector', '7'='NS', '8' = 'NK', '9' = 'NKT',  '10' = 'Treg', '11'='NS', '12'='NS', '13'='NS')
TC.sfj$subcell.id <- Idents(TC.sfj)

# 2. Subset T-cell from NL and Tu tissue
TC_Tu.sfj <- subset(TC.sfj, cells = Tu_TC.bc)
TC_NL.sfj <- subset(TC.sfj, cells = NL_TC.bc)

### Part B. Trajectory
# 1. Subset cell of interes: ex) Tumor CD4+T cells
TC_Tu_CD4.sfj <- subset(TC_Tu.sfj, idents = c("naive CD4", "CD4_memory", "Treg"))
DefaultAssay (TC_Tu_CD4.sfj) <- 'RNA'
pd <- new('AnnotatedDataFrame', data = TC_Tu_CD4.sfj@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(TC_Tu_CD4.sfj), row.names = row.names(TC_Tu_CD4.sfj)))
cds <- newCellDataSet(as(TC_Tu_CD4.sfj@assays$RNA@data, "matrix"), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
cds$clusters <- TC_Tu_CD4.sfj@meta.data$subcell.id 

# 2. Estimate size factors and dispersions
cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)

# 3. Filtering low-quality cells 
set.seed(123) 
cds <- detectGenes(cds, min_expr = 0.1)

# 4. Choosing genes that define progress !
##### Select genes use Seurat::FindAllMarkers function
markers <- FindAllMarkers(object = TC_Tu_CD4.sfj, min.pct = 0.25, thresh.use = 0.25)  
markers <- subset(markers, p_val_adj < 0.05)
nrow(markers)
order.genes <- unique(as.character(markers$gene))

##### Constructing Single Cell Trajectories
cds1 <- setOrderingFilter(cds, order.genes)
cds2 <- reduceDimension(cds = cds1, max_components = 3, method = 'DDRTree') cds2 <- orderCells(cds2) 

##### Plot trajectory
# 1.
png(filename = paste0(Trj.outdir, "CD4_Tu_plot_cell_trajectory.pseudotime.png"), width = 800, height = 800)
plot_cell_trajectory(cds2, color_by = "Pseudotime") 
dev.off()

# 2.
png(filename = paste0(Trj.outdir, "CD4_Tu_plot_cell_trajectory.clusters.png"), width = 800, height = 800)
plot_cell_trajectory(cds2, color_by = "clusters") 
dev.off()

# 3.
png(filename = paste0(Trj.outdir, "CD4_Tu_plot_pseudotime_heatmap.png"), width = 500, height = 1500)
plot_pseudotime_heatmap(cds2[order.genes,], num_clusters = 3, cores = 1, show_rownames = T)
dev.off() 

# 4.
genes <- row.names(subset(fData(cds2), gene_short_name %in% c("TIGIT", "ICOS", "DUSP4", "GNLY", "GZMK")))
cds.sfjbset <- cds2[genes,] 
png(filename = paste0(Trj.outdir, "CD8_plot_genes_in_pseudotime.png", sep=""), width = 800, height = 800)
plot_genes_in_pseudotime(cds.sfjbset, color_by = "clusters")
dev.off()

### Part C. Save results
save.image(paste0(Trj.outdir, "TC_tragectory.RData"))



