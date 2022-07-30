#Set working directory 
setwd("/home/user/Documents/scRNA/datasets")
#Combining two subsets
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(CellChat)
#Read 10x genomics file
#Subset1
sham <- Read10X(data.dir = "sham_adult/")
sham <- CreateSeuratObject(counts = sham, min.cells = 3, min.features  = 200, project = "sham", assay = "RNA")
sham <- NormalizeData(object = sham, normalization.method = "LogNormalize", scale.factor = 10000)
sham <- FindVariableFeatures(object = sham, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = sham))
sham <- ScaleData(object = sham)
#PCA
sham <- RunPCA(object = sham,  npcs = 30, verbose = FALSE)
#Cell clustering
sham <- FindNeighbors(sham, reduction = "pca", dims = 1:20)
sham <- FindClusters(sham, resolution = 0.5, algorithm = 1)
#Creating subset of cdh cells
cdh_sham <- subset(sham, Pecam1 > 1)
cdh_sham <- NormalizeData(object = cdh_sham, normalization.method = "LogNormalize", scale.factor = 10000)
cdh_sham <- FindVariableFeatures(object = cdh_sham, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = cdh_sham))
cdh_sham <- ScaleData(object = cdh_sham)
#PCA
cdh_sham <- RunPCA(object = cdh_sham,  npcs = 30, verbose = FALSE)
cdh_sham <- FindNeighbors(cdh_sham, reduction = "pca", dims = 1:20)
cdh_sham <- FindClusters(cdh_sham, resolution = 0.5, algorithm = 1)

#Subset2
#Read 10x genomics file
mi<- Read10X(data.dir = "mi_adult/")
mi <- CreateSeuratObject(counts = mi, min.cells = 3, min.features  = 200, project = "mi", assay = "RNA")
mi <- NormalizeData(object = mi, normalization.method = "LogNormalize", scale.factor = 10000)
mi <- FindVariableFeatures(object = mi, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = mi))
mi <- ScaleData(object = mi)
#PCA
mi <- RunPCA(object = mi,  npcs = 30, verbose = FALSE)
#Cell clustering
mi <- FindNeighbors(mi, reduction = "pca", dims = 1:20)
mi <- FindClusters(mi, resolution = 0.5, algorithm = 1)
#Creating subset of cdh cells
cdh_mi <- subset(mi, Pecam1 > 1)
cdh_mi <- NormalizeData(object = cdh_mi, normalization.method = "LogNormalize", scale.factor = 10000)
cdh_mi <- FindVariableFeatures(object = cdh_mi, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = cdh_mi))
cdh_mi <- ScaleData(object = cdh_mi)
#PCA
cdh_mi <- RunPCA(object = cdh_mi,  npcs = 30, verbose = FALSE)
cdh_mi <- FindNeighbors(cdh_mi, reduction = "pca", dims = 1:20)
cdh_mi <- FindClusters(cdh_mi, resolution = 0.5, algorithm = 1)

#----------------------###CombineObjects####-------------------------------------------#

hspc.combined <- merge(x = cdh_sham , y = cdh_mi)
hspc.combined <- SplitObject(hspc.combined, split.by = "orig.ident")
hspc.combined <- lapply(X = hspc.combined, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
})
hspc.anchors <- FindIntegrationAnchors(object.list = hspc.combined,dims = 1:15)
hspc.combined <- IntegrateData(anchorset = hspc.anchors, dims = 1:15)

DefaultAssay(hspc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
hspc.combined <- ScaleData(hspc.combined, verbose = FALSE)
hspc.combined <- RunPCA(hspc.combined, npcs = 30, verbose = FALSE)
# UMAP and Clustering
hspc.combined <- RunUMAP(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- FindNeighbors(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- FindClusters(hspc.combined, resolution = 0.5)
hspc.combined <- RunUMAP(hspc.combined, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)

####pecam1 grt1 clusters
hspc.combined[["old.ident"]] <- Idents(object = hspc.combined)
hspc.combined <- RenameIdents(object = hspc.combined, `0` = "aEC1", `1` = "aEC2", `2` = "aEC4", `3` = "aEC5", `4` = "aEC7", `5` = "aEC3", `6` = "Cap EC", `7` = "aEC6",`8` = "SMC", `9` = "EndoMT", `10` = "Cyc aEC")
my_levels <- levels(hspc.combined) <- c("aEC1", "aEC2","aEC3","aEC4","aEC5","aEC6","aEC7","Cyc aEC","Cap EC","EndoMT","SMC")
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
DimPlot(hspc.combined,reduction = "umap",pt.size = 1, cols = c("aEC1" = "greenyellow", "aEC2" = "limegreen", "aEC3" = "seagreen", "aEC4" = "mediumseagreen", "aEC5" = "forestgreen","aEC6" = "yellowgreen","aEC7" = "olivedrab","Cap EC"="orange","Cyc aEC"="red","SMC"="lightgray","EndoMT"="gray")) +NoLegend()
my_levels <- levels(hspc.combined) <- c("SMC","EndoMT","Cap EC","Cyc aEC","aEC7", "aEC6","aEC5","aEC4","aEC3","aEC2","aEC1")
DotPlot(hspc.combined, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Cxcr4", "Gja5", "Jag2", "Bmx","Top2a","Birc5","Cenpa","Mki67","Tyms","Mcm4","Mcm6","Mcm7","Aplnr","Car4","Ssh2","Fmo2","Aqp1","Twist1","Pdgfra","Sox9","Serpine1","Loxl1","Pdgfrb","Notch3","Myh11","Acta2","Rgs5"), assay = "RNA", cols = c("yellow","black")) +RotatedAxis()
