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
cdh_sham <- subset(sham, pecam1 > 1)
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
##subsetting based on orig-.ident
sham <- subset(hspc.combined, subset = orig.ident == "sham")
mi <- subset(hspc.combined, subset = orig.ident == "mi")
#subsetting with required clusters
DimPlot(hspc.combined, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
##dimplots_v2
DimPlot(hspc.combined, reduction = "umap",pt.size = 1.5, cols = c(`0` = "mediumseagreen", `1` = "limegreen", `2` = "seagreen", `3` = "greenyellow", `4` = "forestgreen", `5` = "olivedrab", `6` = "orange", `7` = "yellowgreen",`8` = "darkgray", `9` = "gray", `10`= "orangered",`11`="dimGray"))
DimPlot(sham, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
DimPlot(mi, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
my_levels <- levels(hspc.combined) <- c(3, 0, 5, 2,4, 7, 1, 10, 6, 8, 9, 11)
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
StackedVlnPlot(hspc.combined, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Top2a","Birc5","Cenpa","Aplnr","Car4","Aqp1","Twist1","Pdgfra","Sox9","Pdgfrb","Notch3","Myh11"), color.use = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
##rename clusters
hspc.combined[["old.ident"]] <- Idents(object = hspc.combined)
hspc.combined <- RenameIdents(object = hspc.combined, `3` = "aEC1", `0` = "aEC2", `6` = "Cap EC", `2` = "aEC4", `5` = "aEC3", `4` = "aEC5", `10` = "Cyc aEC", `7` = "aEC6",`8` = "EndoMT", `9` = "SMC1",`11`="SMC2",`1`="aEC7")
my_levels <- levels(hspc.combined) <- c("aEC1", "aEC2","aEC3","aEC4","aEC5","aEC6","aEC7","Cyc aEC","Cap EC","EndoMT","SMC1","SMC2")
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
StackedVlnPlot(hspc.combined, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Top2a","Birc5","Cenpa","Aplnr","Car4","Aqp1","Twist1","Pdgfra","Sox9","Pdgfrb","Notch3","Myh11"), color.use = c(`aEC2` = "cornflowerblue", `aEC7` = "chocolate", `aEC4` = "blue", `aEC1` = "lightskyblue", `aEC5` = "BurlyWood", `aEC#` = "dodgerblue", `Cap EC` = "orange", `aEC6` = "sandybrown",`EndoMT` = "darkgray", `SMC1` = "gray", `Cyc aEC`= "red",`SMC2`="dimGray"))

##Dotplots
my_levels <- levels(hspc.combined) <- c("EndoMT","SMC2","SMC1","Cyc aEC","Cap EC","aEC7", "aEC6","aEC5","aEC4","aEC3","aEC2","aEC1")
DotPlot(hspc.combined, c("tdTomato","Pecam1", "Cldn5", "Fabp4", "Notch1", "Efnb2", "Gja4", "Jag1", "Bmx", "Gja5", "Cxcr4","Sepp1","Aqp1","Aplnr","Car4", "Top2a","Mki67","Birc5","Cenpa","Tyms","Mcm6","Mcm4","Mcm7","Acta2","Tagln","Rgs5","Gm13889","Twist1","Pdgfra","Sox9","Mfap4"),assay = "RNA", cols = c("yellow","black"))
##artery cells subclustered
art<-subset(hspc.combined, idents = c(0,1,2,3,4,5,6,7,10))
art_sham <- subset(art, subset = orig.ident == "sham")
art_mi <- subset(art, subset = orig.ident == "mi")
DimPlot(art, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
DimPlot(art_sham, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
DimPlot(art_mi, reduction = "umap",pt.size = 1.5, cols = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))
my_levels <- levels(art) <- c(3, 0, 5, 2,4, 7, 1, 10, 6)
Idents(art) <- factor(Idents(art), levels= my_levels)
StackedVlnPlot(art, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Top2a","Birc5","Cenpa","Aplnr","Car4","Aqp1"), color.use = c(`0` = "cornflowerblue", `1` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "dodgerblue", `6` = "orange", `7` = "sandybrown",`8` = "darkgray", `9` = "gray", `10`= "red",`11`="dimGray"))

##reclustering
art <- ScaleData(art, verbose = FALSE)
art <- RunPCA(art, npcs = 30, verbose = FALSE)
# UMAP and Clustering
art <- RunUMAP(art, reduction = "pca", dims = 1:20)
art <- FindNeighbors(art, reduction = "pca", dims = 1:20)
art <- FindClusters(art, resolution = 0.5)
art <- RunUMAP(art, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)
DimPlot(art,reduction = "umap",pt.size = 1.5, cols = c(`1` = "dodgerblue", `0` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "sandybrown", `6` = "cornflowerblue", `7` = "peru",`8` = "orange", `9` = "red"))
my_levels <- levels(art) <- c(3, 6, 1,2,4,5,7,0,9,8)
Idents(art) <- factor(Idents(art), levels= my_levels)
StackedVlnPlot(art, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Jag2", "Bmx", "Gja5", "Cxcr4","Top2a","Birc5","Cenpa","Aplnr","Car4","Aqp1"), color.use = c(`1` = "dodgerblue", `0` = "chocolate", `2` = "blue", `3` = "lightskyblue", `4` = "BurlyWood", `5` = "sandybrown", `6` = "cornflowerblue", `7` = "peru",`8` = "orange", `9` = "red"))

no_cells <- table(cluster, cluster$integrated_snn_res.0.5)
no_cells <- write.csv(no_cells, file = "/home/user/Documents/scRNA/plots_cdh5/adult_fig1/no_cells.csv")
#find markers for every cluster compared to all remaining cells
cluster_markers <- FindAllMarkers(object = cluster, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
cluster_markers <- write.csv(cluster_markers, file ="/home/user/Documents/scRNA/plots_cdh5/adult_fig1/markers.csv" )
#Observing slected genes.
FeaturePlot(object = hspc.combined, features = c("Cdh5"), min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1)
FeaturePlot(object = cluster, features = c("tdTomato"), min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1)
FeaturePlot(object = cluster, features = c("Gja5"), min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1)
FeaturePlot(object = cluster, features = c("Aplnr"), min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1)
FeaturePlot(object = cluster, features = c("Prox1"), min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1)
FeaturePlot(object = cluster, features = c("Lyve1"), min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1)
####pecam1 grt1 clusters
hspc.combined[["old.ident"]] <- Idents(object = hspc.combined)
hspc.combined <- RenameIdents(object = hspc.combined, `0` = "aEC1", `1` = "aEC2", `2` = "aEC4", `3` = "aEC5", `4` = "aEC7", `5` = "aEC3", `6` = "Cap EC", `7` = "aEC6",`8` = "SMC", `9` = "EndoMT", `10` = "Cyc aEC")
my_levels <- levels(hspc.combined) <- c("aEC1", "aEC2","aEC3","aEC4","aEC5","aEC6","aEC7","Cyc aEC","Cap EC","EndoMT","SMC")
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
DimPlot(hspc.combined,reduction = "umap",pt.size = 1, cols = c("aEC1"="powderblue", "aEC2"= "lightskyblue","aEC3"="dodgerblue","aEC4" = "burlywood","aEC5"="sandybrown","aEC6" = "peru","aEC7"="olivedrab","Cap EC"="orange","Cyc aEC"="red","SMC"="lightgray","EndoMT"="gray")) +NoLegend()
StackedVlnPlot(hspc.combined, c("tdTomato", "Pecam1", "Cldn5", "Fabp4", "Dll4", "Notch1", "Efnb2", "Hey1", "Gja4", "Cxcr4", "Gja5", "Jag2", "Bmx","Top2a","Birc5","Cenpa","Mki67","Mcm4","Mcm6","Mcm7","Aplnr","Car4","Ssh2","Fmo2","Aqp1","Twist1","Pdgfra","Sox9","Serpine1","Loxl1","Pdgfrb","Notch3","Myh11","Acta2","Rgs5"), color.use = c("aEC1"="powderblue", "aEC2"= "lightskyblue","aEC3"="dodgerblue","aEC4" = "burlywood","aEC5"="sandybrown","aEC6" = "peru","aEC7"="olivedrab","Cap EC"="orange","Cyc aEC"="red","SMC"="lightgray","EndoMT"="gray"))

------------------------------------------#CellCycleSorting of the seurat object-----------------------------------------------
# Basic function to convert human to mouse gene names
library(biomaRt)
library(Seurat)
library(dplyr)
library(cowplot)

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
#imporing data
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
#Running cellcyclesoring
ccs <- CellCycleScoring(object = hspc.combined, g2m.features = m.g2m.genes, s.features = m.s.genes, ident = TRUE)
head(ccs[[]])
ccs <- RunPCA(ccs, features = c(m.s.genes, m.g2m.genes))
DimPlot(ccs, reduction = "pca", group.by = "Phase")
ccs_neo_sham <- write.csv (ccs$Phase, file = "/home/user/Documents/scRNA/plots_cdh5/adult_fig1/ccs.csv")
ccs_neo_sham <- write.csv (ccs$orig.ident, file = "/home/user/Documents/scRNA/plots_cdh5/adult_fig1/ccs_orig_ident.csv")
G1 <- subset(ccs,subset = Phase == "G1")
G2M <- subset(ccs,subset = Phase == "G2M")
S <- subset(ccs,subset = Phase == "S")
DimPlot(G1, pt.size = 1.5, group.by = "orig.ident")
DimPlot(G2M, pt.size = 1.5,group.by = "orig.ident")
DimPlot(S, pt.size = 1.5,group.by = "orig.ident")
--------------------------------#########Trajectory- Importing seurat object to cds format##########---------------------------------------
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(hspc.combined)
cds <- cluster_cells(cds)
#cds <- as.cell_data_set(hspc.combined)
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
#cds <- as.cell_data_set(sham)
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
#cds <- as.cell_data_set(mi)
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
