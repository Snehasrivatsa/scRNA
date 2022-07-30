seurat_clusters <- sham$seurat_clusters
orig_ident<-sham$orig.ident
cyto_bar <- results$CytoTRACE
write.csv(cyto_bar, "/home/user/Documents/scRNA/10x/results_cytotrace_sham.csv")
write.csv(seurat_clusters, "/home/user/Documents/scRNA/10x/seurat_clusters_sham.csv")
write.csv(orig_ident, "/home/user/Documents/scRNA/10x/seurat_idents.csv")
my_data <- read.csv(file.choose())
##neonate
set.seed(1234)
dplyr::sample_n(my_data)
my_data$seurat_clusters <- ordered(my_data$seurat_clusters,
levels = c(0,1,3,4,2,6))
library(dplyr)
group_by(my_data, seurat_clusters, orig_ident) %>%
summarise(
count = n(),
mean = mean(cytotrace_scores
, na.rm = TRUE),
sd = sd(cytotrace_scores, na.rm = TRUE)
)
library("ggpubr")
ggboxplot(my_data, x = c("seurat_clusters","orig_ident"), y = "cytotrace_scores", 
          color = "seurat_clusters", palette = c("YellowGreen", "OliveDrab","darkGreen","forestGreen","Orange","OrangeRed","Gainsboro","Darkgray","DimGray","PowderBlue"),
          order = c(),
          ylab = "cytotrace_score", xlab = "seurat_clusters",
          fill  = "seurat_clusters", notch = T)

##adult
set.seed(1234)
dplyr::sample_n(my_data, 10)
my_data$seurat_clusters <- ordered(my_data$seurat_clusters,
                                   levels = c(3,0,5,2,4,7,1,10,6,8,9,11))
library(dplyr)
group_by(my_data, seurat_clusters) %>%
  summarise(
    count = n(),
    mean = mean(cytotrace_scores
                , na.rm = TRUE),
    sd = sd(cytotrace_scores, na.rm = TRUE)
  )
library("ggpubr")
ggboxplot(my_data, x = "seurat_clusters", y = "cytotrace_scores", 
          color = "seurat_clusters", palette = c("greenyellow", "mediumseagreen","olivedrab","seagreen","forestgreen", "yellowgreen", "limegreen","orangered","orange","darkgray","gray","dimgray"),
          order = c(3,0,5,2,4,7,1,10,6,8,9,11),
          ylab = "cytotrace_score", xlab = "seurat_clusters",
          fill  = "seurat_clusters", notch = T)