# Loading libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(devtools)
library(parallel)

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

# Previewing data
## Data preprocessing
# Two plots to see whether we should normalize data or not.
VlnPlot(brain, features = c("nCount_Spatial", "percent.mt"), ncol = 2)
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Normalizing data with SCTransform
# Very slow, using multisession
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- f1_brain
str(brain)

## Gene expression visualization
# Overlaying molecular data over tissue photo.
sorted_counts <- sort(rowSums(brain_counts[ncrna_feats,]))
split_counts <- split(sorted_counts, ceiling(seq_along(sorted_counts)/4))
SpatialFeaturePlot(brain, features = names(split_counts[[1]]))
SpatialFeaturePlot(brain, features = names(split_counts[[2]]))
SpatialFeaturePlot(brain, features = names(split_counts[[3]]))
SpatialFeaturePlot(brain, features = names(split_counts[[4]]))
SpatialFeaturePlot(brain, features = names(split_counts[[5]]))
SpatialFeaturePlot(brain, features = names(split_counts[[6]]))
SpatialFeaturePlot(brain, features = names(split_counts[[7]]))
SpatialFeaturePlot(brain, features = names(split_counts[[8]]))
SpatialFeaturePlot(brain, features = names(split_counts[[9]]))
SpatialFeaturePlot(brain, features = names(split_counts[[10]]))

# Making plots.
p1 <- SpatialFeaturePlot(brain, features = "Meg3", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Meg3", alpha = c(0.1, 1))
p1 + p2
saveRDS(brain, "pre_reduction.rds")

## Dimensionality reduction, clustering, and visualization
# Dimensionality reduction
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# Visualizing the results of clustering
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
# Singling out clusters
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3)

## Interactive plots
ISpatialDimPlot(brain)

?ISpatialFeaturePlot
ISpatialFeaturePlot(brain, feature = "Meg3")

LinkedDimPlot(brain)

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
saveRDS(des, "de_markers.rds")

SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)


brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "markvariogram")
saveRDS(brain, "find_spatially_variable_features_brain.rds")

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
saveRDS(cortex, "cortex1.rds")

# now remove additional cells, use SpatialDimPlots to visualize what to remove
SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)


p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

saveRDS(cortex, "cortex2.rds")
