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
BiocManager::install("Rsamtools")
install.packages("Signac")
library(Signac)
?Bioconductor
# Set random seed
set.seed(888)
# Loading data
brain <- LoadData("stxBrain", type = "anterior1")

# Finding ncRNA features
# Omitting empty
brain_counts <- brain@assays$Spatial@counts
feats <- rownames(brain_counts[rowSums(brain_counts) > 0,])
ncrnas <- tolower(read.table("./ncRNAs.txt", header=TRUE, sep="	", fill=TRUE)[,2])
ncrnas <- paste0(toupper(substr(ncrnas, 1, 1)), substr(ncrnas, 2, nchar(ncrnas)))
ncrna_feats <- intersect(feats, ncrnas)
ncrna_feats
brain[["percent.mt"]] <- PercentageFeatureSet(brain, features = ncrna_feats)
#brain <- subset(brain, features = ncrna_feats)
str(brain@assays$Spatial@counts@Dimna)

# Previewing data
## Data preprocessing
# Two plots to see whether we should normalize data or not.
VlnPlot(brain, features = c("nCount_Spatial", "percent.mt"), ncol = 2) %>%
    ggsave(file="vlnplot1.svg", plot=.)
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2) %>%
    ggsave(file="nc_vlnplot2.svg", plot=.)

# Normalizing data with SCTransform
# Very slow, using multisession
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
str(brain)

## Gene expression visualization
# Overlaying molecular data over tissue photo.
sorted_counts <- sort(rowSums(brain_counts[ncrna_feats,]))
split_counts <- split(sorted_counts, ceiling(seq_along(sorted_counts)/4))
SpatialFeaturePlot(brain, features = names(split_counts[[1]])) %>%
    ggsave(file="ncRNA1.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[2]])) %>%
    ggsave(file="ncRNA2.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[3]])) %>%
    ggsave(file="ncRNA3.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[4]])) %>%
    ggsave(file="ncRNA4.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[5]])) %>%
    ggsave(file="ncRNA5.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[6]])) %>%
    ggsave(file="ncRNA6.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[7]])) %>%
    ggsave(file="ncRNA7.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[8]])) %>%
    ggsave(file="ncRNA8.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[9]])) %>%
    ggsave(file="ncRNA9.svg", plot=.)

SpatialFeaturePlot(brain, features = names(split_counts[[10]])) %>%
    ggsave(file="ncRNA10.svg", plot=.)


# Making plots.
p1 <- SpatialFeaturePlot(brain, features = "Meg3", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Meg3", alpha = c(0.2, 1))
p1 + p2 %>%
    ggsave(file="Meg3.svg", plot=.)
saveRDS(brain, "pre_reduction.rds")
brain <- readRDS("pre_reduction.rds")
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
