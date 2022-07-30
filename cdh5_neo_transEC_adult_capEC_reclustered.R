--------------------------##NEONATE##---------------------------------------------
##subsetting cap
cap <- subset(hspc.combined, idents =c("trans-EC"))
# Run the standard workflow for visualization and clustering
cap <- ScaleData(cap, verbose = FALSE)
cap <- RunPCA(cap, npcs = 30, verbose = FALSE)
# UMAP and Clustering
cap <- RunUMAP(cap, reduction = "pca", dims = 1:20)
cap <- FindNeighbors(cap, reduction = "pca", dims = 1:20)
cap <- FindClusters(cap, resolution = 0.5)
cap <- RunUMAP(cap, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)
##rename clusters
cap[["old.ident"]] <- Idents(object = cap)
cap <- RenameIdents(object = cap, `0` = "artery-like", `1` = "capEC", `2` = "venous-cap", `3` = "EndoMT-like")
##Dimplot
my_cols <- c("artery-like"="#FAC898", "capEC"="#FF5F1F","venous-cap"="#F28F1C","EndoMT-like"="#B87333")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
DimPlot(cap,cols = my_cols2)
##subset based on orig.ident
mi<-subset(cap, orig.ident == "mi")
sham<-subset(cap, orig.ident == "sham")
##Dotplot
my_levels <- levels(cap) <- c("EndoMT-like","venous-cap","capEC","artery-like")
DotPlot(cap, c("tdTomato","Gja5","Cxcr4","Gja4","Dll4","Hey1","Aplnr","Gpihbp1","Pcdh17","Lamb1","Emcn","Vwf","Nr2f2","Plvap","Myl2","Tnnt2","Tnni3","Myl3"), assay = "RNA", cols = c("yellow","black"))


--------------------------------------------##ADULT##------------------------------------------
##subsetting cap
cap <- subset(hspc.combined, idents =c("Cap EC"))
# Run the standard workflow for visualization and clustering
cap <- ScaleData(cap, verbose = FALSE)
cap <- RunPCA(cap, npcs = 30, verbose = FALSE)
# UMAP and Clustering
cap <- RunUMAP(cap, reduction = "pca", dims = 1:20)
cap <- FindNeighbors(cap, reduction = "pca", dims = 1:20)
cap <- FindClusters(cap, resolution = 0.5)
cap <- RunUMAP(cap, n.neighbors = 10, dims = 1:20, spread = 2, min.dist = 0.3)
##rename clusters
cap[["old.ident"]] <- Idents(object = cap)
cap <- RenameIdents(object = cap, `0` = "capEC", `1` = "art-capEC", `2` = "venous-capEC")
##Dimplot
my_cols <- c("art-capEC"="#FAC898", "capEC"="#FF5F1F","venous-capEC"="#F28F1C")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
DimPlot(cap,cols = my_cols2)
##subset based on orig.ident
mi<-subset(cap, orig.ident == "mi")
sham<-subset(cap, orig.ident == "sham")
##Dotplot
my_levels <- levels(cap) <- c("venous-capEC","capEC","art-capEC")
DotPlot(cap, c("tdTomato","Cxcr4","Gja5","Gja4","Dll4","Aplnr","Car4","Aqp1","Rbp7","Emcn","Vwf","Nr2f2","Plvap","Top2a","Birc5","Mki67","Cenpa","Tyms","Mcm4","Mcm6","Mcm7"), assay = "RNA", cols = c("yellow","black"))
