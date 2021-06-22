### Part A. set up
library(Seurat)
Assign.outdir = "./Results/"
ITG.suj <- readRDS("{  /   }") # get previous dataset
ITG.list <- SplitObject(ITG.suj, split.by = "orig.ident")

# Define function
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_logFC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_logFC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

# Hiracheal clustering and ClusterTree
DimPlot(ITG.suj, reduction = "umap", label = T)
PlotClusterTree(BuildClusterTree(ITG.suj))


### Part B. Assign cluster
# 1. Keep previous labels
ITG.suj$subcluster.id <- Idents(ITG.suj)

# 2. Tissue.id
# Tu tissue origin cell barcodes
Tu.bc <- c(colnames(ITG.list$Case_01_Tu),colnames(ITG.list$Case_02_Tu), colnames(ITG.list$Case_03_Tu), colnames(ITG.list$Case_04_Tu), colnames(ITG.list$Case_05_Tu), colnames(ITG.list$Case_06_Tu))
# NL tissue origin cell barcodes
NL.bc <- c(colnames(ITG.list$Case_01_NL),colnames(ITG.list$Case_02_NL), colnames(ITG.list$Case_03_NL), colnames(ITG.list$Case_04_NL), colnames(ITG.list$Case_05_NL))
## Assign
Idents(ITG.suj, cells = NL.bc) <- "NL"
Idents(ITG.suj, cells = Tu.bc) <- "Tu"
ITG.suj$tissue.id <- Idents(ITG.suj)
# Reset subcluster.id
Idents(ITG.suj) <- ITG.suj$subcluster.id


### Part C. Get Conserved Markers
DefaultAssay(ITG.suj) <- 'RNA'
for (i in levels(ITG.suj)) {
  FCM_i <- FindConservedMarkers(ITG.suj, ident.1 = i, grouping.var = 'orig.ident', assay = 'RNA')
  write.table(FCM_i, file = paste0(ITG_MK.outdir, "ITG_", i, "_cnv.mk"), quote = F, row.names = TRUE, sep = "\t")
}


### Part D. FindAllmarkers
ITG_All.mk <- FindAllMarkers(ITG.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(ITG_All.mk, file = paste0(ITG_MK.outdir, "/Integ_ITG_All.mk"), quote = F, row.names = FALSE, sep = "\t")

ITG_All.pw <- FindAllMarkers(ITG.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
ITG_All <- merge(ITG_All.mk, ITG_All.pw, by="row.names", all.x = F, all.y = F)
ITG_All <- subset(ITG_All, select = c('Row.names', 'cluster.x', 'avg_logFC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
ITG_All <- ITG_All[rev(order(ITG_All$cluster.x, ITG_All$avg_logFC.x, ITG_All$power)), ]
write.table(ITG_All, file = paste0(ITG_MK.outdir, "/Integ_ITG_All.mp"), quote = F, sep = "\t")

### Part E. Find each cluster markers. 
for (idx in levels(ITG.suj)) {
  a <- FindMarkers(ITG.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE)
  b <- FindMarkers(ITG.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("ITG_",idx,".mp"), c)
  write.table(c, file = paste0(ITG_MK.outdir, "ITG_", idx, ".mp"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}


### Part F. assign celltype id #################################################
ITG.sfj <- ITG.suj

# Cancer (CA) barcode
# Case_01
Case_01_CA_1.bc <- readLines('/home/yschang/CODES/Case_01_CA_1.bc')
Case_01_CA_2.bc <- readLines('/home/yschang/CODES/Case_01_CA_2.bc')

# Case_02
Case_02_CA_1.bc <- readLines('/home/yschang/CODES/Case_02_CA_1.bc')
Case_02_CA_2.bc <- readLines('/home/yschang/CODES/Case_02_CA_2.bc')

# Case_03
Case_03_CA_1.bc <- readLines('/home/yschang/CODES/Case_03_CA_1.bc')

# Case_04
Case_04_CA_1.bc <- readLines('/home/yschang/CODES/Case_04_CA_1.bc')
Case_04_CA_2.bc <- readLines('/home/yschang/CODES/Case_04_CA_2.bc')

# Case_05
Case_05_CA_1.bc <- readLines('/home/yschang/CODES/Case_05_CA_1.bc')
Case_05_CA_2.bc <- readLines('/home/yschang/CODES/Case_05_CA_2.bc')

# Case_06
Case_06_CA_1.bc <- readLines('/home/yschang/CODES/Case_06_CA_1.bc')
Case_06_CA_2.bc <- readLines('/home/yschang/CODES/Case_06_CA_2.bc')


# Cancer cell barcodes
CA_1.bc<- c(Case_01_CA_1.bc, Case_02_CA_1.bc, Case_03_CA_1.bc, Case_04_CA_1.bc, Case_05_CA_1.bc, Case_06_CA_1.bc)
CA_2.bc<- c(Case_01_CA_2.bc, Case_02_CA_2.bc, Case_04_CA_2.bc, Case_05_CA_2.bc)
CA.bc <- c(CA_1.bc, CA_2.bc)


# CA(Neoplastic)
ITG.sfj <- SetIdent(ITG.sfj, cells = CA.bc, value = "CA")

# EP
EP_temp.bc <- WhichCells(subset(ITG.sfj, idents = c('2', '3', '6', '10', '12', '13', '17')))
EP.bc <- setdiff(EP_temp.bc, CA.bc)
ITG.sfj <- SetIdent(ITG.sfj, cells = EP.bc, value = "EP")

# TC
ITG.sfj <- RenameIdents(ITG.sfj, '0' = 'TC', '1' = 'TC')

# BC
ITG.sfj <- RenameIdents(ITG.sfj, '14' = 'BC')

# MY
ITG.sfj <- RenameIdents(ITG.sfj, '4' = 'MY', '5' = 'MY', '7' = 'MY', '15' = 'MY', '16' ='MY')

# MA
ITG.sfj <- RenameIdents(ITG.sfj, '8' = 'MA')

# EC
ITG.sfj <- RenameIdents(ITG.sfj, '9' = 'EC')

# FB
ITG.sfj <- RenameIdents(ITG.sfj, '11' = 'FB')


# Assign idents and visualize
ITG.sfj$celltype.id <- Idents(ITG.sfj)
DimPlot(ITG.sfj, reduction = "umap", label = T, pt.size = 1.0)

# Save
save.image(paste0(Assign.outdir, "{ITG_}.RData"))
gc(); q()
