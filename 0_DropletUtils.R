# Set-Up
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(Seurat)

# Get dataset and crete seurat object
dir_Exp1 <- "/media/yschang/V/MySamples/outs/count/filtered_feature_bc_matrix/"
agg_raw.data <- Read10X(dir_Exp1)
agg.sbj <- CreateSeuratObject(counts = agg_raw.data, project = "agg")

# Seurat object to SCE
agg.sce <- as.SingleCellExperiment(agg.sbj)

# Manipulate data
rownames(agg.sce) <- uniquifyFeatureNames(rowData(agg.sce)$ID, rowData(agg.sce)$Symbol)
colnames(agg.sce) <- agg.sce$Barcode
my.counts <- counts(agg.sce)
br.out <- barcodeRanks(my.counts)

# generate plots
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))

# Estimating empty droplets
set.seed(100)
e.out <- emptyDrops(my.counts)
# In conclusion, out data do hot have empty droplets

# Removing droplets using FDR value
is.cell <- e.out$FDR <= 0.05 # <= 0.01
sum(is.cell, na.rm=TRUE)

is.cell[is.na(is.cell)] <- FALSE

names(is.cell) <- colnames(agg.sce)
agg.sce$cells_kept <- is.cell

agg.sce <- agg.sce[, agg.sce$cells_kept == T]
agg.sce

### In conclusion, out data do hot have empty droplets