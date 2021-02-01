#QC, Nomralization, plots using seurat for scRNA data.
 
#clean workspace
rm(list=ls())
 
# Load required packages
library(Seurat)
library(dplyr)
library(cowplot)
 
#Set working directory 
setwd("/home/soumyashreelab/Data")
 
#Read 10x genomics file
mice_neonate <- Read10X(data.dir = "mi_adult/")
 
#Create Seurat object
mi_neo <- CreateSeuratObject(counts = mice_neonate, min.cells = 3, min.features  = 200, project = "10X_mi_neo", assay = "RNA")

#Data Normalization
mi_neo <- NormalizeData(object = mi_neo, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable genes across scRNA
mi_neo <- FindVariableFeatures(object = mi_neo, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
#To view the file with the variable genes:
##head(x = HVFInfo(object = mi_neo))

#Data Scaling
mi_neo <- ScaleData(object = mi_neo)


#Linear dimensional reduction:
#PCA
mi_neo <- RunPCA(object = mi_neo,  npcs = 30, verbose = FALSE)
#Cell clustering
mi_neo <- FindNeighbors(mi_neo, reduction = "pca", dims = 1:20)

#Computing nearest neighbourhood graph
mi_neo <- FindClusters(mi_neo, resolution = 0.5, algorithm = 1)

#QC remove mictochondrial RNA
mito.genes <- grep(pattern = "^MT-", x = rownames(mi_neo@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(mi_neo@assays[["RNA"]][mito.genes, ])/Matrix::colSums(mi_neo@assays[["RNA"]])
mi_neo <- AddMetaData(object = mi_neo, metadata = percent.mito, col.name = "percent.mito")  
mi_neo$percent.mito <- percent.mito
VlnPlot(object = mi_neo, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

#Scatter plot of mitochondrial RNA
par(mfrow = c(1, 2))
FeatureScatter(object = mi_neo, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = mi_neo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
mi_neo <- subset(x = mi_neo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito >  -Inf & percent.mito < 0.05)


#Scatter plot
DimHeatmap(object = mi_neo, reduction = "pca", cells = 200, balanced = TRUE)

#Plot PCA
DimPlot(object = mi_neo, reduction = "pca")
#Jack straw plot
#mi_neo <- ScoreJackStraw(object = mi_neo, dims = 1:20, reduction = "pca")
#JackStrawPlot(object = mi_neo, dims = 1:20, reduction = "pca")

#Elbow plot.
#ElbowPlot(object = mi_neo)


#Non-linear dimensional reduction (tSNE)
mi_neo <- RunTSNE(object = mi_neo, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = mi_neo, reduction = "tsne")

#UMAP
mi_neo <- RunUMAP(mi_neo, reduction = "pca", dims = 1:20)
DimPlot(mi_neo, reduction = "umap", label = TRUE)


#Finding deferentially expressed genes
#find markers for every cluster compared to all remaining cells, report only the positive ones
mi_neo.markers <- FindAllMarkers(object = mi_neo, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
mi_neo.markers <- write.csv(mi_neo.markers, file ="/home/soumyashreelab/Data/markers_mi_adult.csv" )

#Featured plot of a endothelial markers
FeaturePlot(object = mi_neo, features = c("Cxcl12","Kdr", "EMCN", "Pecam1","Jag1","Gja5","Cdh5","tdTomato", label= TRUE))

#feature plot cardiomyocytes markers
FeaturePlot(object = mi_neo, features = c("TNNT2","MYOM","Nkx2-5","Tbx5","Isl1", label= TRUE))

#feature plot smooth muscle markers  
FeaturePlot(object = mi_neo, features = c("Acta2","Myh11","PDGFRB","Notch3", label= TRUE))

#feature plot macrophage markers 
FeaturePlot(object = mi_neo, features = c("Cd11","Cd14","Cd86","CX3CR1","Mafb","Mhc2","Spp1","Arg1","Chi3l3","TGFB1", label= TRUE))


#Visualize protein levels on RNA clusters
#Adding expression data to either the counts, data, 
#or scale.data slots can be done with SetAssayData. 
#New data must have the same cells in the same order 
#as the current expression data. Data added to counts or data
#must have the same features as the current expression data.
#Normalize with the assay data.
mi_neo.adt <- GetAssayData(object = mi_neo, slot = 'data')
mi_neo[["ADT"]] <- CreateAssayObject(counts = mi_neo.adt)
mi_neo <- NormalizeData(mi_neo, assay = "ADT", normalization.method = "CLR")
mi_neo <- ScaleData(mi_neo, assay = "ADT")
FeaturePlot(mi_neo, features = "tdTomato", min.cutoff = "q05", max.cutoff = "q95", ncol = 4)


#cells per cluster
cells_per_cluster <- table(Idents(mi_neo))
cells_per_cluster <- as.data.frame(cells_per_cluster) %>%
rename(Cluter = Var1,
       Number_cells = Freq)
cells_per_cluster <- write.csv(cells_per_cluster, file = "/home/soumyashreelab/Data/cells_per_cluster_mi_adult.csv")
