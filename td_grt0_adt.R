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
cdh_sham <- subset(sham, tdTomato > 0)
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
cdh_mi <- subset(mi, tdTomato > 0)
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
hspc.combined <- RenameIdents(object = hspc.combined, `0` = "aEC1", `1` = "aEC2", `6` = "aEC3", `3` = "aEC4", `2` = "aEC5", `5` = "aEC6", `8` = "Cyc aEC",`7`="SMC",`4` = "EndoMT")
my_levels <- levels(hspc.combined) <- c("aEC1","aEC2","aEC3","aEC4","aEC5","aEC6","Cyc aEC","SMC","EndoMT")
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
DimPlot(hspc.combined,pt.size = 1.5, cols = c("aEC1"="lightskyblue","aEC2"="Dodgerblue","aEC3"="blue","aEC4"="burlywood","aEC5"="sandybrown","aEC6"="chocolate","Cyc aEC"="red","SMC"="gray","EndoMT"="dimgray"))
StackedVlnPlot(hspc.combined, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Top2a","Birc5","Cenpa","Mki67","Mcm7","Myh11","Notch3","Cspg4","Myocd","Pdgfrb","Serpine1","Cdh2","Mfap4","Twist1","Sox9"), color.use = c("aEC1"="lightskyblue","aEC2"="Dodgerblue","aEC3"="blue","aEC4"="burlywood","aEC5"="sandybrown","aEC6"="chocolate","Cyc aEC"="red","SMC"="gray","EndoMT"="dimgray"))

##subsetting and reclustering based on orig.ident
sham <- subset(hspc.combined, subset = orig.ident == "sham")
mi <- subset(hspc.combined, subset = orig.ident == "mi")

sham <- ScaleData(sham, verbose = FALSE)
sham <- RunPCA(sham, npcs = 30, verbose = FALSE)
# UMAP and Clustering
sham <- RunUMAP(sham, reduction = "pca", dims = 1:20)
sham <- FindNeighbors(sham, reduction = "pca", dims = 1:20)
sham <- FindClusters(sham, resolution = 0.5)
sham <- RunUMAP(sham, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)

sham <- RenameIdents(object = sham, `10` = "aEC1", `6` = "aEC2", `0` = "aEC3", `3` = "aEC4", `1` = "aEC5", `5` = "aEC6",`2`="aEC7",`4`="aEC8", `9` = "Cyc aEC",`8`="SMC",`7` = "EndoMT")
my_levels <- levels(sham) <- c("aEC1","aEC2","aEC3","aEC4","aEC5","aEC6","aEC7","aEC8","Cyc aEC","SMC","EndoMT")
Idents(sham) <- factor(Idents(sham), levels= my_levels)
DimPlot(sham,pt.size = 1, cols = c("aEC1"="powderblue","aEC2"="lightskyblue","aEC3"="deepskyblue","aEC4"="yellowgreen","aEC5"="forestgreen","aEC6"="burlywood","aEC7"="sandybrown","aEC8"="chocolate","Cyc aEC"="red","SMC"="gray","EndoMT"="dimgray"))
##Dimplot v2:
DimPlot(sham,pt.size = 1, cols = c("aEC1"="mediumseagreen","aEC2"="palegreen","aEC3"="seagreen","aEC4"="yellowgreen","aEC5"="limegreen","aEC6"="olivedrab","aEC7"="greenyellow","aEC8"="forestgreen","Cyc aEC"="red","SMC"="gray","EndoMT"="dimgray"))
StackedVlnPlot(sham, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Top2a","Mki67","Birc5", "Tyms","Mcm6","Mcm4",'Mcm7',"Myh11","Notch3","Cspg4","Speg","Pdgfrb","Col1a1","Cdh2","Mfap4","Twist1","Sox9"), color.use = c("aEC1"="powderblue","aEC2"="lightskyblue","aEC3"="deepskyblue","aEC4"="yellowgreen","aEC5"="forestgreen","aEC6"="burlywood","aEC7"="sandybrown","aEC8"="chocolate","Cyc aEC"="red","SMC"="gray","EndoMT"="dimgray"))

##Dotplot
my_levels <- levels(sham) <- c("EndoMT","SMC","Cyc aEC","aEC8","aEC7","aEC6","aEC5","aEC4","aEC3","aEC2","aEC1")
DotPlot(sham, c("tdTomato","Cdh5","Pecam1","Tek","Cldn5","Gja5","Cxcr4","Bmx","Gja4","Notch1","Efnb2","Jag1","Top2a","Mki67","Birc5", "Tyms","Mcm6","Mcm4",'Mcm7',"Myh11","Notch3","Cspg4","Pdgfrb","Col1a1","Mfap4","Twist1","Sox9","Smad3","Zeb2","Snai1","Snai2","Gata6","Loxl2","Malat1","Acvrl1","Tgfbr1"),assay = "RNA",cols = c("yellow","black"))

mi <- ScaleData(mi, verbose = FALSE)
mi <- RunPCA(mi, npcs = 30, verbose = FALSE)
# UMAP and Clustering
mi <- RunUMAP(mi, reduction = "pca", dims = 1:20)
mi <- FindNeighbors(mi, reduction = "pca", dims = 1:20)
mi <- FindClusters(mi, resolution = 0.5)
mi <- RunUMAP(mi, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)

mi <- RenameIdents(object = mi, `5` = "aEC1", `0` = "aEC2", `3` = "aEC3", `2` = "aEC4", `6` = "EndoMT1",`7`="SMC",`1`="EndoMT2",`4` = "EndoMT3")
my_levels <- levels(mi) <- c("aEC1","aEC2","aEC3","aEC4","SMC","EndoMT1","EndoMT2","EndoMT3")
Idents(mi) <- factor(Idents(mi), levels= my_levels)
DimPlot(mi,pt.size = 1.5, cols = c("aEC1"="burlywood","aEC2"="sandybrown","aEC3"="peru","aEC4"="forestgreen","SMC"="gainsboro","EndoMT1"="gray","EndoMT2"="darkgray","EndoMT3"="dimgray"))
##Dimplot v2
DimPlot(mi,pt.size = 1.5, cols = c("aEC1"="olivedrab","aEC2"="greenyellow","aEC3"="limegreen","aEC4"="forestgreen","SMC"="gainsboro","EndoMT1"="gray","EndoMT2"="darkgray","EndoMT3"="dimgray"))
StackedVlnPlot(mi, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Myh11","Notch3","Cspg4","Myocd","Cnn1","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7',"Loxl1","Cdh2","Mfap4","Twist1","Sox9"), color.use = c("aEC1"="burlywood","aEC2"="sandybrown","aEC3"="peru","aEC4"="forestgreen","SMC"="gainsboro","EndoMT1"="gray","EndoMT2"="darkgray","EndoMT3"="dimgray"))

##DotPlot
my_levels <- levels(mi) <- c("EndoMT3","EndoMT2","EndoMT1","SMC","aEC4","aEC3","aEC2","aEC1")
DotPlot(mi, c("tdTomato","Cdh5","Pecam1","Tek","Cldn5","Gja5","Cxcr4","Bmx","Gja4","Notch1","Efnb2","Jag1","Pdgfrb","Myh11","Notch3","Cspg4","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7',"Pdgfra","Mfap4","Twist1","Sox9","Smad3","Zeb2","Snai1","Snai2","Gata6","Loxl2","Malat1","Acvrl1","Tgfbr1"), assay = "RNA", cols = c("yellow","black"))
--------------------------------#########Trajectory- Importing seurat object to cds format##########---------------------------------------
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(hspc.combined)
cds <- cluster_cells(cds)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_color = "green",
           trajectory_graph_segment_size = 2)
##sham
cds <- as.cell_data_set(sham)
cds <- cluster_cells(cds)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_color = "green",
           trajectory_graph_segment_size = 2)
##mi
cds <- as.cell_data_set(mi)
cds <- cluster_cells(cds)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_color = "green",
           trajectory_graph_segment_size = 2)
---------------#############CytoTRACE##############------------------
library(CytoTRACE)
results <- hspc.combined[["RNA"]]@counts
results <- as.data.frame(results)
umap <- as.data.frame(Embeddings(hspc.combined, reduction = "umap"))
results <- subset(hspc.combined, subset())
results <- CytoTRACE(results)
plotCytoTRACE(results, emb = umap)
##sham
results <- sham[["RNA"]]@counts
results <- as.data.frame(results)
umap <- as.data.frame(Embeddings(sham, reduction = "umap"))
results <- CytoTRACE(results)
plotCytoTRACE(results, emb = umap)
##mi
results <- mi[["RNA"]]@counts
results <- as.data.frame(results)
umap <- as.data.frame(Embeddings(mi, reduction = "umap"))
results <- CytoTRACE(results)
plotCytoTRACE(results, emb = umap)
