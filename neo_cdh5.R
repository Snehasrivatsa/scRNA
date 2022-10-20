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
sham <- Read10X(data.dir = "sham_neonate/")
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
cdh_sham <- subset(sham, Cdh5 > 1)
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
cdh_mi <- subset(mi, Cdh5 > 1)
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

#QC remove mictochondrial RNA
mito.genes <- grep(pattern = "^mt-", x = rownames(hspc.combined@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(hspc.combined@assays[["RNA"]][mito.genes, ])/Matrix::colSums(hspc.combined@assays[["RNA"]])
hspc.combined <- AddMetaData(object = hspc.combined, metadata = percent.mito, col.name = "percent.mito")  
hspc.combined$percent.mito <- percent.mito
VlnPlot(object = hspc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

# UMAP and Clustering
hspc.combined <- RunUMAP(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- FindNeighbors(hspc.combined, reduction = "pca", dims = 1:20)
hspc.combined <- FindClusters(hspc.combined, resolution = 0.5)
hspc.combined <- RunUMAP(hspc.combined, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)
##rename clusters
hspc.combined[["old.ident"]] <- Idents(object = hspc.combined)
hspc.combined <- RenameIdents(object = hspc.combined, `0` = "aEC1", `1` = "aEC2", `2` = "trans-EC", `3` = "aEC3", `4` = "aEC4", `5` = "SMC2", `6` = "cyc-aEC", `7` = "endo-EC",`8` = "endoMT", `9` = "SMC1")
my_levels <- levels(hspc.combined) <- c("aEC1", "aEC2","aEC3","aEC4","trans-EC","cyc-aEC","endoMT","SMC1","SMC2","endo-EC")
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
##subclustering  based on orig.idents
sham <- subset(hspc.combined, subset = orig.ident == "sham")
mi <- subset(hspc.combined, subset = orig.ident == "mi")
##td cells
td = GetAssayData(object = hspc.combined, 
                               assay = "RNA", slot = "data")["tdTomato",]
pos_ids = names(which(td>0))
td = subset(hspc.combined,cells=pos_ids)
td_sham <- subset(td, subset = orig.ident == "sham")
td_mi <- subset(td, subset = orig.ident == "mi")
##artery cells subclustered
art<-subset(hspc.combined, idents = c("aEC1","aEC2","aEC3","aEC4","trans-EC","cyc-aEC"))
art_sham <- subset(art, subset = orig.ident == "sham")
art_mi <- subset(art, subset = orig.ident == "mi")
##reclustering
art <- ScaleData(art, verbose = FALSE)
art <- RunPCA(art, npcs = 30, verbose = FALSE)
# UMAP and Clustering
art <- RunUMAP(art, reduction = "pca", dims = 1:20)
art <- FindNeighbors(art, reduction = "pca", dims = 1:20)
art <- FindClusters(art, resolution = 0.5)
art <- RunUMAP(art, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)
##colour-coordinating
DimPlot(art, pt.size = 1.5, cols = c("4"="yellowgreen","3"="olivedrab","0"="forestgreen","2"="darkgreen","6"="springgreen","1"="orange","7"="orangered","5"="darkred"))
art[["old.ident"]] <- Idents(object = art)
my_levels <- levels(art) <- c(4,3,0,2,6,1,7,5)
Idents(art) <- factor(Idents(art), levels= my_levels)
StackedVlnPlot(art, features = c("tdTomato","Hey1","Dll4","Notch1","Efnb2","Pecam1",  "Cldn5","Fabp4","Sox17","Car4","Aqp7","Fabp3","Pltp","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7'), color.use = c("4"="yellowgreen","3"="olivedrab","0"="forestgreen","2"="darkgreen","6"="springgreen","1"="orange","7"="orangered","5"="darkred"))

DimPlot(hspc.combined, reduction = "umap", group.by = "orig.ident", pt.size = 1)
DimPlot(hspc.combined, reduction = "umap", pt.size = 1.5, cols = c(`0`= "YellowGreen", `1`="OliveDrab",`3`="darkGreen",`4`="forestGreen",`2`="Orange",`6`="OrangeRed",`8`="Gainsboro",`9`="Darkgray",`5`="DimGray",`7`="PowderBlue"))
DimPlot(sham, reduction = "umap", pt.size = 1.5,cols = c(`0`= "YellowGreen", `1`="OliveDrab",`3`="DarkGreen",`4`="ForestGreen",`2`="Orange",`6`="OrangeRed",`8`="Gainsboro",`9`="Darkgray",`5`="DimGray",`7`="PowderBlue"))
DimPlot(mi, reduction = "umap", pt.size = 1.5,label = T,cols = c(`0`= "YellowGreen", `1`="OliveDrab",`3`="DarkGreen",`4`="ForestGreen",`2`="Orange",`6`="OrangeRed",`8`="Gainsboro",`9`="Darkgray",`5`="DimGray",`7`="PowderBlue"))
no_cells <-  table(hspc.combined@active.ident, hspc.combined@meta.data$orig.ident)
no_cells <- write.csv(no_cells, file = "/home/user/Documents/scRNA/Neo_fig1/no_cells.csv")
#find markers for every cluster compared to all remaining cells
cluster_markers <- FindAllMarkers(object = cluster, min.pct = 0.25, thresh.use = 0.25)
cluster_markers <- write.csv(cluster_markers, file ="/home/user/Documents/scRNA/plots_cdh5/adult_3,11_removed/markers.csv" )
#Cells per cluster
cells_per_cluster <- table(Idents(hspc.combined))
cells_per_cluster <- as.data.frame(cells_per_cluster) %>%
  rename(Cluster = Var1,
         Number_cells = Freq)
cells_per_cluster <- write.csv(cells_per_cluster, file = "/home/user/Documents/scRNA/plots_cdh5/adult_3,11_removed/CellsPerCluster_clustered.csv")
#Observing slected genes.
FeaturePlot(object = hspc.combined, features = c("Cdh5","tdTomato","Gja5","Aplnr"),ncol=2, min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1.5) 
FeaturePlot(object = hspc.combined, features = c("Prox1","Lyve1","Prox1","Lyve1"),ncol=2, min.cutoff = 0, max.cutoff = 3, cols = c("white","purple"), pt.size = 1.5)
StackedVlnPlot(hspc.combined, features = c("tdTomato","Pecam1","Cldn5","Fabp4","Hey1","Dll4","Notch1","Efnb2","Pecam1",  "Cldn5","Fabp4","Sox17","Car4","Aqp7","Fabp3","Pltp","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7',"Sox9","Pdgfra","Mfap4","Twist1","Myh11","Pdgfrb","Notch3","Myocd","Npr3","Tgfbi","Apoc1","Adgrg6","Cdh11"), color.use = c("aEC1"= "YellowGreen", "aEC2"="OliveDrab","aEC3"="ForestGreen","aEC4"="DarkGreen","cap-EC"="Orange","cyc-aEC"="OrangeRed","endoMT"="Gainsboro","SMC1"="Darkgray","SMC2"="DimGray","endo-EC"="PowderBlue"))
##Dot plot
my_levels <- levels(hspc.combined) <- c("endo-EC","endoMT","SMC2","SMC1","cyc-aEC","trans-EC","aEC4", "aEC3","aEC2","aEC1")
DotPlot(hspc.combined, features = c("tdTomato","Pecam1","Cldn5","Fabp4", "Notch1", "Efnb2", "Gja4", "Jag1", "Bmx", "Gja5", "Cxcr4", "Rgcc","Car4","Aqp7","Fabp3","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7',"Myh11","Pdgfrb","Notch3","Myocd","Sox9","Pdgfra","Mfap4","Twist1","Npr3","Cdh11","Tgfbi","Adgrg6"), cols = c("yellow","black"), assay = "RNA")

##Pecam1_grt1 colours
DimPlot(hspc.combined,cols = c("0"= "Yellowgreen", "3"="olivedrab","5"= "forestgreen","1"="darkgreen","8"="springgreen","2"="orange","11"="brown", "7"="orangered","4"="Dimgray","10"="darkgray","9"="gray","6"= "powderblue","12"="deepskyblue"),pt.size = 1.5)
##rename clusters
hspc.combined[["old.ident"]] <- Idents(object = hspc.combined)
hspc.combined <- RenameIdents(object = hspc.combined, `0` = "aEC1", `3` = "aEC2", `5` = "aEC3", `1` = "aEC4", `8` = "aEC5", `2` = "trans-EC1", `11` = "trans-EC2", `7` = "cyc-aEC",`9` = "endoMT", `10` = "SMC1", `4` = "SMC2", `6` = "EndoEC1", `12` = "EndoEC2")
my_levels <- levels(hspc.combined) <- c("aEC1", "aEC2","aEC3","aEC4","cap-EC","cyc-aEC","endoMT","SMC1","SMC2","endo-EC")
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
my_levels <- levels(hspc.combined) <- c("EndoEC2","EndoEC1", "SMC2","SMC1","endoMT","cyc-aEC","trans-EC2","trans-EC1","aEC5","aEC4","aEC3","aEC2","aEC1")
 ##pecam1,stackedvln
DotPlot(hspc.combined, features = c("tdTomato","Pecam1",  "Cldn5","Fabp4","Sox17","Hey1","Dll4","Notch1","Efnb2","Car4","Aqp7","Fabp3","Pltp","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7',"Sox9","Pdgfra","Mfap4","Twist1","Myh11","Pdgfrb","Notch3","Myocd", "Npr3","Cdh11","Tgfbi", "Adgrg6"),assay = "RNA",cols=c("yellow","black")) + RotatedAxis()
##pecam1 clusters
hspc.combined[["old.ident"]] <- Idents(object = hspc.combined)
my_levels <- levels(hspc.combined) <- c(0,3,5,1,8,2,11,7,9,10,4,6,12)
Idents(hspc.combined) <- factor(Idents(hspc.combined), levels= my_levels)
##pecam1,stackedvln
StackedVlnPlot(hspc.combined, features = c("tdTomato","Hey1","Dll4","Notch1","Efnb2","Pecam1",  "Cldn5","Fabp4","Sox17","Car4","Aqp7","Fabp3","Pltp","Top2a","Mki67","Birc5","Cenpa", "Tyms","Mcm6","Mcm4",'Mcm7',"Sox9","Pdgfra","Mfap4","Twist1","Myh11","Pdgfrb","Notch3","Myocd", "Npr3","Cdh11","Tgfbi", "Adgrg6"), color.use = c("0"= "Yellowgreen", "3"="olivedrab","5"= "forestgreen","1"="darkgreen","8"="springgreen","2"="orange","11"="brown", "7"="orangered","4"="Dimgray","10"="darkgray","9"="gray","6"= "powderblue","12"="deepskyblue"))

--------------------###########CellCycleSorting of the seurat object######################--------------------
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
#DimPlot(ccs, reduction = "pca", group.by = "Phase")
ccs_neo_sham <- write.csv (ccs$Phase, file = "/home/user/Documents/scRNA/Neo_fig1/phase.csv")
ccs_neo_sham <- write.csv (ccs$orig.ident, file = "/home/user/Documents/scRNA/Neo_fig1/ccs.csv")
G1 <- subset(ccs,subset = Phase == "G1")
G2M <- subset(ccs,subset = Phase == "G2M")
S <- subset(ccs,subset = Phase == "S")
DimPlot(G1, pt.size = 1.5, group.by = "orig.ident")
DimPlot(G2M, pt.size = 1.5,group.by = "orig.ident")
DimPlot(S, pt.size = 1.5,group.by = "orig.ident")
#DimPlot(ccs, reduction = "umap",split.by = "Phase",group.by = "orig.ident")
--------------------------------#########Trajectory- Importing seurat object to cds format##########---------------------------------------
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(mi)
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = TRUE)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = TRUE)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)
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

