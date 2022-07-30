#Set working directory 
setwd("/home/user/suraj")
#Combining two subsets
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(SeuratDisk)
##sham
Convert("adata_sham.h5ad", dest = "h5seurat", overwrite = TRUE)
sham <- LoadH5Seurat("adata_sham.h5seurat")
sham[["orig.ident"]] <- "sham"
sham <- NormalizeData(object = sham, normalization.method = "LogNormalize", scale.factor = 10000)
sham <- FindVariableFeatures(object = sham, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = sham))
# Run the standard workflow for visualization and clustering
sham <- ScaleData(sham, verbose = FALSE)
sham <- RunPCA(object = sham,  npcs = 30, verbose = FALSE)
#Cell clustering
sham <- FindNeighbors(sham, reduction = "pca", dims = 1:20)
sham <- FindClusters(sham, resolution = 0.5, algorithm = 1)
#Creating subset of cdh cells
cdh_sham <- subset(sham, Cdh5 > 1)
cdh_sham <- NormalizeData(object = cdh_sham, normalization.method = "LogNormalize", scale.factor = 10000)
cdh_sham <- FindVariableFeatures(object = cdh_sham, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = cdh_sham))
cdh_sham <- ScaleData(object = cdh_sham)
#PCA
cdh_sham <- RunPCA(object = cdh_sham,  npcs = 30, verbose = FALSE)
cdh_sham <- FindNeighbors(cdh_sham, reduction = "pca", dims = 1:20)
cdh_sham <- FindClusters(cdh_sham, resolution = 0.5, algorithm = 1)


##mi
Convert("adata_MI.h5ad", dest = "h5seurat", overwrite = TRUE)
mi <- LoadH5Seurat("adata_MI.h5seurat")
mi[["orig.ident"]] <- "mi"
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
cdh_mi <- subset(mi, Cdh5 > 1)
cdh_mi <- NormalizeData(object = cdh_mi, normalization.method = "LogNormalize", scale.factor = 10000)
cdh_mi <- FindVariableFeatures(object = cdh_mi, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = cdh_mi))
cdh_mi <- ScaleData(object = cdh_mi)
#PCA
cdh_mi <- RunPCA(object = cdh_mi,  npcs = 30, verbose = FALSE)
cdh_mi <- FindNeighbors(cdh_mi, reduction = "pca", dims = 1:20)
cdh_mi <- FindClusters(cdh_mi, resolution = 0.5, algorithm = 1)

#----------------------###CombineObjects####-------------------------------------------##
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
##subsetting based on orig-.ident
sham <- subset(hspc.combined, subset = orig.ident == "sham")
mi <- subset(hspc.combined, subset = orig.ident == "mi")
#subsetting with required clusters
DimPlot(hspc.combined, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
##dimplots_v2
DimPlot(hspc.combined, reduction = "umap",pt.size = 1.5, cols = c(`0` = "mediumseagreen", `1` = "limegreen", `2` = "seagreen", `3` = "greenyellow", `4` = "forestgreen", `5` = "olivedrab", `6` = "orange", `7` = "yellowgreen",`8` = "darkgray", `9` = "gray", `10`= "orangered",`11`="dimGray"))
DimPlot(sham, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
DimPlot(mi, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))

#strip
hspc.combined_strip <- subset(hspc.combined, idents = c(0,1,2,3,4,7))
mi_strip <- subset(mi, idents = c(0,1,2,3,4,7))
sham_strip <- subset(sham, idents = c(0,1,2,3,4,7))


#RNA velocity

library(velocyto.R)
#install SeuratWrappers
install.packages("devtools")
library(devtools)
devtools::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

#merged_strip velocity
s_cellranger_orig <- RunVelocity(object = hspc.combined_strip, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- c(`0`= "YellowGreen", `1`="Orange",`2`="OliveDrab",`3`="forestGreen",`4`="darkgreen",`7`="OrangeRed")
names(x = ident.colors) <- levels(x = s_cellranger_orig)
cell.colors <- ident.colors[Idents(object = s_cellranger_orig)]
names(x = cell.colors) <- colnames(x = s_cellranger_orig)

show.velocity.on.embedding.cor(emb = Embeddings(object = s_cellranger_orig, reduction = "umap"), vel = Tool(object = s_cellranger_orig, 
                                                                                                            slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

#mi_strip velocity

s_cellranger_orig1 <- RunVelocity(object = mi_strip, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- c(`0`= "YellowGreen", `1`="Orange",`2`="OliveDrab",`3`="forestGreen",`4`="darkgreen",`7`="OrangeRed")
names(x = ident.colors) <- levels(x = s_cellranger_orig1)
cell.colors <- ident.colors[Idents(object = s_cellranger_orig1)]
names(x = cell.colors) <- colnames(x = s_cellranger_orig1)

show.velocity.on.embedding.cor(emb = Embeddings(object = s_cellranger_orig1, reduction = "umap"), vel = Tool(object = s_cellranger_orig1, 
                                                                                                            slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

#sham_strip velocity

s_cellranger_orig2 <- RunVelocity(object = sham_strip, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- c(`0`= "YellowGreen", `1`="Orange",`2`="OliveDrab",`3`="forestGreen",`4`="darkgreen",`7`="OrangeRed")
names(x = ident.colors) <- levels(x = s_cellranger_orig2)
cell.colors <- ident.colors[Idents(object = s_cellranger_orig2)]
names(x = cell.colors) <- colnames(x = s_cellranger_orig2)

show.velocity.on.embedding.cor(emb = Embeddings(object = s_cellranger_orig2, reduction = "umap"), vel = Tool(object = s_cellranger_orig2, 
                                                                                                            slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)



