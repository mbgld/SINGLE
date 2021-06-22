# set up
library(Seurat)
ITG.outdir = "{  }/results/"

################################################################################
### COMBINING DATASETS #########################################################
################################################################################
### Part A. Call datasets
# Case_01 dataset
# NL
Case_01_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/Case_01/A_data_NL/filtered_feature_bc_matrix/"); Case_01_NL.sbj <- CreateSeuratObject(counts = Case_01_NL_raw.data, project = "Case_01_NL", min.cells = 3, min.features = 200); rm(Case_01_NL_raw.data)
# Tu
Case_01_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/Case_01/A_data_Tu/filtered_feature_bc_matrix/"); Case_01_Tu.sbj <- CreateSeuratObject(counts = Case_01_Tu_raw.data, project = "Case_01_Tu", min.cells = 3, min.features = 200); rm(Case_01_Tu_raw.data)

# Case_02 dataset
# NL
Case_02_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/Case_02/A_data_NL/filtered_feature_bc_matrix/"); Case_02_NL.sbj <- CreateSeuratObject(counts = Case_02_NL_raw.data, project = "Case_02_NL", min.cells = 3, min.features = 200); rm(Case_02_NL_raw.data)
# Tu
Case_02_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/Case_02/A_data_Tu/filtered_feature_bc_matrix/"); Case_02_Tu.sbj <- CreateSeuratObject(counts = Case_02_Tu_raw.data, project = "Case_02_Tu", min.cells = 3, min.features = 200); rm(Case_02_Tu_raw.data)

# Case_03 dataset
#	NL
Case_03_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/Case_03/A_data_NL/filtered_feature_bc_matrix/"); Case_03_NL.sbj <- CreateSeuratObject(counts = Case_03_NL_raw.data, project = "Case_03_NL", min.cells = 3, min.features = 200); rm(Case_03_NL_raw.data)
# Tu
Case_03_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/Case_03/A_data_Tu/filtered_feature_bc_matrix/"); Case_03_Tu.sbj <- CreateSeuratObject(counts = Case_03_Tu_raw.data, project = "Case_03_Tu", min.cells = 3, min.features = 200); rm(Case_03_Tu_raw.data)

# Case_04 dataset
# NL
Case_04_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/Case_04/A_data_NL/filtered_feature_bc_matrix/"); Case_04_NL.sbj <- CreateSeuratObject(counts = Case_04_NL_raw.data, project = "Case_04_NL", min.cells = 3, min.features = 200); rm(Case_04_NL_raw.data)
# Tu
Case_04_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/Case_04/A_data_Tu/filtered_feature_bc_matrix/"); Case_04_Tu.sbj <- CreateSeuratObject(counts = Case_04_Tu_raw.data, project = "Case_04_Tu", min.cells = 3, min.features = 200); rm(Case_04_Tu_raw.data)

# Case_05 dataset
# NL
Case_05_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/Case_05/A_data_NL/filtered_feature_bc_matrix/"); Case_05_NL.sbj <- CreateSeuratObject(counts = Case_05_NL_raw.data, project = "Case_05_NL", min.cells = 3, min.features = 200); rm(Case_05_NL_raw.data)
# TU
Case_05_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/Case_05/A_data_Tu/filtered_feature_bc_matrix/"); Case_05_Tu.sbj <- CreateSeuratObject(counts = Case_05_Tu_raw.data, project = "Case_05_Tu", min.cells = 3, min.features = 200); rm(Case_05_Tu_raw.data)

# Case_06 dataset
# Tu
Case_06_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/Case_06/A_data_Tu/filtered_feature_bc_matrix/"); Case_06_Tu.sbj <- CreateSeuratObject(counts = Case_06_Tu_raw.data, project = "Case_06_Tu", min.cells = 3, min.features = 200); rm(Case_06_Tu_raw.data)


### Part B. Merge
Combi.sbj <- merge(Case_01_NL.sbj, y = c( Case_01_Tu.sbj, Case_02_NL.sbj, Case_02_Tu.sbj, Case_06_Tu.sbj, Case_04_NL.sbj, Case_04_Tu.sbj, Case_03_NL.sbj, Case_03_Tu.sbj, Case_05_NL.sbj, Case_05_Tu.sbj), add.cell.ids = c('Case_01_NL', 'Case_01_Tu', 'Case_02_NL', 'Case_02_Tu', 'Case_04_NL', 'Case_04_Tu', 'Case_03_NL', 'Case_03_Tu', 'Case_05_NL', 'Case_05_Tu', 'Case_06_Tu')); rm(Case_01_NL.sbj, Case_01_Tu.sbj, Case_02_NL.sbj, Case_02_Tu.sbj, Case_04_NL.sbj, Case_04_Tu.sbj, Case_03_NL.sbj, Case_03_Tu.sbj, Case_05_NL.sbj, Case_05_Tu.sbj, Case_06_Tu.sbj)


### Part C. QC
# Get mito percentage, QC, and visualize
Combi.sbj[["percent.mt"]] <- PercentageFeatureSet(Combi.sbj, pattern = "^MT-")
Combi.sqj <- subset(Combi.sbj, subset = percent.mt < 20 & nCount_RNA >100 & nCount_RNA < 150000 & nFeature_RNA > 200 & nFeature_RNA < 10000)
VlnPlot(Combi.sqj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Combi.sqj, feature1 = "nCount_RNA", feature2 = "percent.mt"); plot2 <- FeatureScatter(Combi.sqj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"); CombinePlots(plots=list(plot1, plot2))


### Part D. LogNormalization and Scaling
Combi.snj = NormalizeData(Combi.sqj, normalization.method = "LogNormalize", scale.factor = 10000)
Combi.snj <- FindVariableFeatures(Combi.snj, selection.method = "vst", nfeatures = 2000); v.genes <- VariableFeatures(Combi.snj)
LabelPoints(plot = VariableFeaturePlot(Combi.snj), points = v.genes[1:15],repel = T)
Combi.ssj <- ScaleData(Combi.snj, features = rownames(Combi.snj))


### Part E. Dimension Reduction
# 1. Linear Dimensional Reduction by PCA. 
Combi.spj <- RunPCA(Combi.ssj, features = v.genes)
print(Combi.spj[["pca"]], dims = 1:5, nfeatures=5) 
VizDimLoadings(Combi.spj, dims = 1:2, reduction ="pca")
DimPlot(Combi.spj, reduction = 'pca')
DimHeatmap(Combi.spj, dims = 1, cells = 500, balanced =TRUE)
DimHeatmap(Combi.spj, dims = 1:15, cells = 500, balanced =TRUE)

# 2. NON-LINEAR DIMENSION REDUCTION 
ElbowPlot(Combi.spj)
Combi.spj <- JackStraw(Combi.spj, num.replicate = 100)
Combi.spj <- ScoreJackStraw(Combi.spj, dims = 1:20)
JackStrawPlot(Combi.spj, dims = 1:20)
Combi.spj <- FindNeighbors(Combi.spj, dims = 1:10)  
Combi.spj <- FindClusters(Combi.spj, resolution = 0.2)
Combi.suj <- RunUMAP(Combi.spj, dims = 1:10) 

# 3. Visualize
DimPlot(Combi.suj, reduction = "umap", label = T)
DimPlot(Combi.suj, group.by= 'orig.ident', label = T)

##### END of merging Every cases ###############################################

################################################################################
### INTEGRATING DATASETS #######################################################
################################################################################


### Part A. Dataset Preprocessing 
dataset.list <- SplitObject(Combi.suj, split.by = 'orig.ident')
dataset.list <- lapply(dataset.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})
Anchors <- FindIntegrationAnchors(object.list = dataset.list, dims=1:30)

### Part B. Integration
ITG.sbj <- IntegrateData(anchorset=Anchors, dims=1:30)
DefaultAssay(ITG.sbj) <- 'integrated'

### Part C. Scaling
ITG.ssj <- ScaleData(ITG.sbj, features = rownames(ITG.sbj))

### Part D. Dimensional Reduction

# 1. Linear Dimensional Reduction by PCA. 
ITG.spj <- RunPCA(ITG.ssj, npcs = 30, verbose = FALSE); print(ITG.spj[["pca"]], dims = 1:5, nfeatures=5)
VizDimLoadings(ITG.spj, dims = 1:2, reduction ="pca"); DimPlot(ITG.spj, reduction = 'pca'); DimHeatmap(ITG.spj, dims = 1, cells = 500, balanced =TRUE); DimHeatmap(ITG.spj, dims = 1:15, cells = 500, balanced =TRUE)

# 2.NON-LINEAR DIMENSION REDUCTION
ElbowPlot(ITG.spj)
ITG.spj <- JackStraw(ITG.spj, num.replicate = 100)
ITG.spj <- ScoreJackStraw(ITG.spj, dims = 1:20)
JackStrawPlot(ITG.spj, dims = 1:20)
ITG.spj <- FindNeighbors(ITG.spj, dims = 1:10) 
ITG.spj <- FindClusters(ITG.spj, resolution = 0.2)
ITG.suj <- RunUMAP(ITG.spj, reduction = 'pca',  dims = 1:10)
DimPlot(ITG.suj, reduction = "umap", label = T)
DimPlot(ITG.suj, group.by= 'orig.ident', label = T)

##################################################
save.image(paste0(ITG.outdir, "{   }.RData"))
gc()
q()
