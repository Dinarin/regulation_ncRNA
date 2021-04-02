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
str(brain)

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
f1_brain <- future(SCTransform(brain, assay = "Spatial", verbose = FALSE))

saveRDS(f1_brain, file = "SCTransform_result.rds")

resolved(f1_brain)
brain <- value(f1_brain)

## Gene expression visualization
# Overlaying molecular data over tissue photo.
SpatialFeaturePlot(brain, features = ncrna_feats)

empty_feats <- c("Nron", "Tsix", "Chaer1", "Plut", "Hottip", "Hotair", "Mir155hg")
very_low_feats <- c("Dnm3os", "Paupar", "Xist", "Rbakdn", "Hotairm1", "Carmn")
low_feats <- c("Mhrt", "Flicr")

SpatialFeaturePlot(brain, features = very_low_feats)


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


brain <- readRDS("brain_before_finding_spatial_features.rds")
des <- markers <- value(readRDS("markers.rds"))

SpatialFeaturePlot(object = brain, features = rownames(des <- markers)[1:3], alpha = c(0.1, 1), ncol = 3)


# Very slow, using future
f3_brain <- future(FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "markvariogram"))
resolved(f3_brain)
brain <- value(f3_brain)

f_top.features <- future(head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6))
top.features <- value(f_top.features)
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
