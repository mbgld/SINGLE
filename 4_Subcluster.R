### Set_up
library(Seurat)
Subcluster.outdir = './Results/'

### Load previous work
# get dataset from previous work: Splitted into each major cell cluster

# Define useful fxs
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_log2FC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

### Part A. Reclustering
DefaultAssay(Cluster.sbj) <- 'RNA'
Cluster.sbj <- FindVariableFeatures(Cluster.sbj, selection.method = "vst", nfeatures = 2000, assay = 'RNA'); v.genes <- VariableFeatures(Cluster.sbj)
DefaultAssay(Cluster.sbj) <- 'integrated'
Cluster.ssj <- ScaleData(object = Cluster.sbj, features = rownames(Cluster.sbj)); rm(Cluster.sbj)
Cluster.spj <- RunPCA(object = Cluster.ssj, features = v.genes)

# Estimate dimension
ElbowPlot(object = Cluster.spj)
Cluster.spj <- JackStraw(Cluster.spj, num.replicate = 100); Cluster.spj <- ScoreJackStraw(Cluster.spj, dims = 1:20); JackStrawPlot(Cluster.spj, dims = 1:20)

# UMAP
Cluster.suj <- RunUMAP(Cluster.spj, dims = 1:5)
Cluster.suj <- FindNeighbors(Cluster.suj, dims = 1:10)
Cluster.suj <- FindClusters(Cluster.suj, resolution = 0.2)
DimPlot(Cluster.suj, reduction = "umap", label = T, group.by = 'orig.ident')
DimPlot(Cluster.suj, reduction = "umap", label = T, pt.size = 1.0)
DimPlot(Cluster.suj, reduction = "umap", label = F, group.by = 'tissue.id', pt.size = 1.0)

### Part B. Find markers
DefaultAssay(Cluster.suj) <- 'RNA'

# 1. FindAllmarkers
Cluster_All.mk <- FindAllMarkers(Cluster.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(Cluster_All.mk, file = paste0(Cluster.outdir, "/Integ_Cluster_All.mk"), quote = F, row.names = FALSE, sep = "\t")
Cluster_All.pw <- FindAllMarkers(Cluster.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
Cluster_All <- merge(Cluster_All.mk, Cluster_All.pw, by="row.names", all.x = F, all.y = F)
Cluster_All <- subset(Cluster_All, select = c('Row.names', 'cluster.x', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
Cluster_All <- Cluster_All[rev(order(Cluster_All$cluster.x, Cluster_All$avg_log2FC.x, Cluster_All$power)), ]
write.table(Cluster_All, file = paste0(Cluster.outdir, "/Integ_Cluster_All.mp"), quote = F, sep = "\t")

# 2. Find each cluster markers. 
for (idx in levels(Cluster.suj)) {
  a <- FindMarkers(Cluster.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE)
  b <- FindMarkers(Cluster.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("Cluster_sc_",idx,".mp"), c)
  write.table(c, file = paste0(Cluster.outdir, "/Integ_Cluster_", idx, ".mp"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}


### Part C. Counts
# 1. count each cluster
for (i in levels(Cluster.suj)){
  print(paste0("Cluster_", i, " length: ", length(WhichCells(Cluster.suj, idents = i))))
  rm(i)
}

# 2. count by tissue.id
for (i in levels(Cluster.suj)){
  print(paste0("Cluster_", i, " in NL Tissue count: ", length(intersect(YS_NL.bc, WhichCells(Cluster.suj, idents = i)))))
  print(paste0("Cluster_", i, " in Tu Tissue count: ", length(intersect(YS_Tu.bc, WhichCells(Cluster.suj, idents = i)))))
  rm(i)
}

# 3. count by project.id
Cluster.list <- SplitObject(Cluster.suj, split.by = "orig.ident")
## KJA
for (i in levels(Cluster.suj)){
  print(paste0("KJA_NL cluster_", i, " length: ", length(intersect(colnames(Cluster.list$KJA_NL), WhichCells(Cluster.suj, idents = i)))))
  print(paste0("KJA_Tu cluster_", i, " length: ", length(intersect(colnames(Cluster.list$KJA_Tumor), WhichCells(Cluster.suj, idents = i)))))
  rm(i)
}


### Part D. Assing idents
# 1. Keep subcluster.id
Cluster.sfj$subcluster.id <- Idents(Cluster.suj)

# 2. subcell.id
Cluster.sfj <- RenameIdents(Cluster.sfj, '0' = 'Anti-inflammatory AM', '1' = 'TAM', '2' = 'CD1c+ DC', '3' = 'Proinflammatory mo-Mac', '4' = 'NS', '5' = 'Proinflammatory mo-Mac', '6' = 'NS', '7' = 'CD163+CD14+DCs', '8'= 'Cycling AM',  '9' = 'Proinflammatory mo-Mac', '10'= 'Activated DC')
Cluster.sfj$subcell.id <- Idents(Cluster.sfj)

### Part  Save
save.image(paste0(Cluster.outdir, "1_Integ_Cluster.RData"))
gc();q()


