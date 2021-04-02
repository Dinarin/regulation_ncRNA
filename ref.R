# Loading libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(devtools)
library(harmony)

# check the current active plan
plan()
# change the current plan to access parallelization
plan(multisession, workers = 8)
availableCores()
plan()
options(future.globals.maxSize = 7 * 1024^3)

# Loading data
brain <- LoadData("stxBrain", type = "anterior1")

# Previewing data
## Data preprocessing
# Two plots to see whether we should normalize data or not.
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Normalizing data with SCTransform
# Very slow, using multicore on Linux
f1_brain <- future(SCTransform(brain, assay = "Spatial", verbose = FALSE))
saveRDS(f1_brain, file = "SCTransform_result.rds")
resolved(f1_brain)
brain <- value(f1_brain)

## Gene expression visualization
# Overlaying molecular data over tissue photo.
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
# Making plots.
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2


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


ISpatialFeaturePlot(brain, features = "Ttr")

LinkedDimPlot(brain)

f_des <- f_markers <- future(FindMarkers(brain, ident.1 = 5, ident.2 = 6))
resolved(f_des)
saveRDS(f_markers, file = "markers.rds")
saveRDS(brain, file = "brain_before_finding_spatial_features.rds")
des <- markers <- value(f_des)
SpatialFeaturePlot(object = brain, features = rownames(des <- markers)[1:3], alpha = c(0.1, 1), ncol = 3)


# Very slow, using multicore
f3_brain <- future(FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "markvariogram"))
resolved(f3_brain)
brain <- value(f3_brain)

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image <- imagerow > 400 |
# image <- imagecol < 150))
cortex <- subset(cortex, anterior1 <- imagerow > 400 | anterior1 <- imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1 <- imagerow > 275 & anterior1 <- imagecol > 370, invert =
		  TRUE)
cortex <- subset(cortex, anterior1 <- imagerow > 250 & anterior1 <- imagecol > 440, invert =
		  TRUE)


p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2










allen <- reference <- readRDS("./data/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells
# this speeds up SCTransform dramatically with no loss in performance

library(dplyr)
allen <- reference <- SCTransform(allen <- reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)


# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen <- reference, group.by = "subclass", label = TRUE)


anchors <- FindTransferAnchors(reference = allen <- reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen <- reference$subclass, prediction.assay = TRUE,
				      weight.reduction = cortex[["pca"]])
cortex[["predictions"]] <- predictions.assay


DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)


cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram",
					    features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)



SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
					    "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)


brain.merge <- merge(brain, brain2)                                                   

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)


DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))


SpatialDimPlot(brain.merge)


SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))


slide.seq <- LoadData("ssHippo")
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log <- nCount <- Spatial <- log(slide.seq$nCount <- Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap <- plots(plot1, plot2)
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c(1, 6, 13)), facet.highlight = TRUE)

