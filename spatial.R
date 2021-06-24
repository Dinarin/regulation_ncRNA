# Loading libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


library(future)
library(devtools)
library(parallel)
library(svglite)
library(Signac)

# Set random seed
set.seed(888)

# Setting parallel
plan("multiprocess", workers = 4)

# Loading data
brain_a1 <- LoadData("stxBrain", type = "anterior1")
brain_p1 <- LoadData("stxBrain", type = "posterior1")

?stxBrain
# Finding ncRNA features
ncrnas <- tolower(read.table("./ncRNAs.txt", header=TRUE, sep="	", fill=TRUE)[,2])
ncrnas <- paste0(toupper(substr(ncrnas, 1, 1)), substr(ncrnas, 2, nchar(ncrnas)))

# Omitting empty
finding_relevant_feats <- function (brain_data, target_feats) {
  brain_counts <- brain_data[["Spatial"]]@counts
  feats <- intersect(rownames(brain_counts), target_feats)
  return(feats)
}

brain_a1_ncrna_feats <- finding_relevant_feats(brain_a1, ncrnas)
brain_p1_ncrna_feats <- finding_relevant_feats(brain_p1, ncrnas)

## Previewing data
# Data preprocessing
# Two plots to see whether we should normalize data or not.

plot1 <- VlnPlot(brain_a1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain_a1, features = "nCount_Spatial")
wrap_plots(plot1, plot2) %>%
  ggsave(file="a1_vlnplot.svg", plot=.)

plot1 <- VlnPlot(brain_p1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain_p1, features = "nCount_Spatial")
wrap_plots(plot1, plot2) %>%
  ggsave(file="p1_vlnplot.svg", plot=.)

## Integrating data
brain_list <- list(anterior1 = brain_a1, posterior1 = brain_p1)

# Running SCT on both datasets
brain_list <- lapply(brain_list, SCTransform, assay = "Spatial", method = "poisson")

# Setting maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # setting allowed size to 2K MiB


brain_features <- SelectIntegrationFeatures(brain_list, nfeatures = 3000, verbose = FALSE)
brain_list <- PrepSCTIntegration(object.list = brain_list, anchor.features = brain_features, verbose = FALSE)
brain_anchors <- FindIntegrationAnchors(object.list = brain_list, normalization.method = "SCT", verbose = FALSE, anchor.features = brain_features)
brain_integrated <- IntegrateData(anchorset = brain_anchors, normalization.method = "SCT", verbose = FALSE)

rm(brain_anchors, brain_list)
gc()

# Overlaying molecular data over tissue photo.
a1_sorted_counts <- sort(rowSums(brain_a1[["Spatial"]]@counts[brain_a1_ncrna_feats,]))
a1_split_counts <- split(sorted_counts, ceiling(seq_along(sorted_counts)/4))
SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[1]])) %>%
    ggsave(file="a1_ncRNA1.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[2]])) %>%
    ggsave(file="a1_ncRNA2.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[3]])) %>%
    ggsave(file="a1_ncRNA3.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[4]])) %>%
    ggsave(file="a1_ncRNA4.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[5]])) %>%
    ggsave(file="a1_ncRNA5.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[6]])) %>%
    ggsave(file="a1_ncRNA6.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[7]])) %>%
    ggsave(file="a1_ncRNA7.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[8]])) %>%
    ggsave(file="a1_ncRNA8.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[9]])) %>%
    ggsave(file="a1_ncRNA9.svg", plot=.)

SpatialFeaturePlot(brain_a1, features = names(a1_split_counts[[10]])) %>%
    ggsave(file="a1_ncRNA10.svg", plot=.)

p1_sorted_counts <- sort(rowSums(brain_p1[["Spatial"]]@counts[brain_p1_ncrna_feats,]))
p1_split_counts <- split(p1_sorted_counts, ceiling(seq_along(p1_sorted_counts)/4))
SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[1]])) %>%
    ggsave(file="p1_ncRNA1.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[2]])) %>%
    ggsave(file="p1_ncRNA2.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[3]])) %>%
    ggsave(file="p1_ncRNA3.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[4]])) %>%
    ggsave(file="p1_ncRNA4.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[5]])) %>%
    ggsave(file="p1_ncRNA5.svg", plot=.)
?SpatialFeaturePlot
SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[6]])) %>%
    ggsave(file="p1_ncRNA6.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[7]])) %>%
    ggsave(file="p1_ncRNA7.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[8]])) %>%
    ggsave(file="p1_ncRNA8.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[9]])) %>%
    ggsave(file="p1_ncRNA9.svg", plot=.)

SpatialFeaturePlot(brain_p1, features = names(p1_split_counts[[10]])) %>%
    ggsave(file="p1_ncRNA10.svg", plot=.)


# Making plots.
p1 <- SpatialFeaturePlot(brain, features = "Meg3", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Meg3", alpha = c(0.2, 1))
p1 + p2 %>%
    ggsave(file="Meg3.svg", plot=.)
load("after_SCT.RData")

## Dimensionality reduction, clustering, and visualization
# Dimensionality reduction
?RunPCA
??svd
brain <- RunSVD(brain, assay = "SCT", features = ncrna_feats, verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# Visualizing the results of clustering
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
wrap_plots(p1 + p2) %>%
    ggsave(file="clusters.svg", plot=., width=12, height=6)
# Singling out clusters
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3) %>%
    ggsave(file="SDplot.svg", plot=.)

## Interactive plots
ISpatialDimPlot(brain)

?ISpatialFeaturePlot
ISpatialFeaturePlot(brain, feature = "Meg3")

LinkedDimPlot(brain)

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
saveRDS(des, "de_markers.rds")
de_markers <- readRDS("markers.rds")

SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3) %>%
    ggsave(file="SFmarkers_plot.svg", plot=.)


brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "markvariogram")
saveRDS(brain, "find_spatially_variable_features_brain.rds")
brain <- readRDS("find_features.rds")

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)

SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1)) %>%
ggsave(file="SFplot.svg", plot=.)


cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
saveRDS(cortex, "cortex1.rds")
cortex <- readRDS("cortex1.rds")

# now remove additional cells, use SpatialDimPlots to visualize what to remove
#SpatialDimPlot(cortex, cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)


p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2 %>%
ggsave(file="cortex_plot.svg", plot=.)

saveRDS(cortex, "cortex2.rds")
