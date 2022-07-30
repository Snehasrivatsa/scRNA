#Set working directory 
setwd("/home/user/Documents/scRNA/datasets")
#Combining two subsets
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
# library(CellChat)
#Read 10x genomics file
#Subset1
sham <- Read10X(data.dir = "sham_neonate/")
sham <- CreateSeuratObject(counts = sham, min.cells = 3, min.features  = 200, project = "sham", assay = "RNA")
sham <- NormalizeData(object = sham, normalization.method = "LogNormalize", scale.factor = 10000)
sham <- FindVariableFeatures(object = sham, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = sham))
sham <- ScaleData(object = sham)
#PCA  ``
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
mi<- Read10X(data.dir = "mi_neonate/")
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
##subclustering  based on orig.idents
sham <- subset(hspc.combined, subset = orig.ident == "sham")
mi <- subset(hspc.combined, subset = orig.ident == "mi")
##markers
DotPlot(hspc.combined, features = c("tdTomato","Aplnr","Cxcr4","Gja5","Gja4","Hey1","Dll4","Notch1","Efnb2","Bmx","Jag1","Jag2","Cdh5","Pecam1","Cldn5","Fabp4","Sox17"), assay = "RNA", cols = c("yellow","black"))
##CMs
DotPlot(hspc.combined, features = c("Nkx2-5","Tbx5","Tnnt2"), assay = "RNA", cols = c("yellow","black"))
##macrophages
DotPlot(hspc.combined, features = c("Spp1","Cx3cr1"), assay = "RNA", cols = c("yellow","black"))
##mes
DotPlot(hspc.combined, features = c("Myl2","Actc1","Tnni3","Mb","Tnnc1","Tnnt2","Myh6","Myl3","Pln","Nrp2","Cd300lg","Mest"), assay = "RNA", cols = c("yellow","black"))
##cap
DotPlot(hspc.combined, features = c("Car4","Aqp1","Adm","Mycn","Cxcl12","Rbp7","Mgll","Ly6c1","Aqp7","Btnl9","Pdgfd","Sema7a","Clec2d","Unc5b"), assay = "RNA", cols = c("yellow","black"))
##cyc
DotPlot(hspc.combined, c("Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7'), assay = "RNA", cols = c("yellow","black"))
##art
DotPlot(hspc.combined, c("Fbln5","Stmn2","Glul","Fos","Id1","Gadd45g","Btg2","Ier2","Klf4","Alpl","Junb","Dusp1","Cyr61","Mecom","Hspa1a","Klf2","Jun","Plk2","Crip1","Slc6a6"), assay = "RNA", cols = c("yellow","black"))
##LEC
DotPlot(hspc.combined, c("Prox1","Lyve1","Gpm6a","Lbp", "Thy1","Mmrn1","Fgl2","Fth1","Nts"), assay = "RNA", cols = c("yellow","black"))
##SMCs
DotPlot(hspc.combined, c("Myh11","Pdgfrb","Notch3","Myocd","Col26a1", "Speg"), assay = "RNA", cols = c("yellow","black"))
##EndoMTs
DotPlot(hspc.combined, c("Serpine1","Loxl1","Col1a1","Tagln","Mfap4","Twist1", "Pdgfra", "Sox9","Cdh2","Col26a1"), assay = "RNA", cols = c("yellow","black"))
##large vein
DotPlot(hspc.combined, features = c("Nr2f2","Vwf","Mgp","Cfh","Apoe","Cpe","Cytl1"), assay = "RNA", cols = c("yellow","black"))
DotPlot(hspc.combined,c("Emcn","Bgn","Plvap","Dcn","Ctsh","Rbp1","Npr3","H19","Tm4sf1","Igfbp4","Fabp5","Tmem108","Id2","Cgnl1","Clu","Hmcn1","Il6st","Ece1","Atp1b3","Tmem176b"), assay = "RNA", cols = c("yellow","black"))
##venous capillary
DotPlot(hspc.combined, features = c("Vcam1","Ier3", "Pltp","Pi16","Fmo2"), assay = "RNA", cols = c("yellow","black"))
##capillary
DotPlot(hspc.combined, features = c("Kdr","Ssh2", "Endou","Fam57b","Kifc3"), assay = "RNA", cols = c("yellow","black"))
##endocardium
DotPlot(hspc.combined, c("Nfatc1","Npr3","Gata5","Tgfbi","Apoc1"), assay = "RNA", cols = c("yellow","black"))
##fibroblasts
DotPlot(hspc.combined, c("Col1a1","Col1a2","Col5a1","Loxl1","Lum","Fbln1","Fbln2","Cd34","Pdgfra"), assay = "RNA", cols = c("yellow","black"))