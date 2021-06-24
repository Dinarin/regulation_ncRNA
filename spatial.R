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

# Set random seed
set.seed(888)

# Setting parallel
plan("multiprocess", workers = 4)

# Loading data
brain_a1 <- LoadData("stxBrain", type = "anterior1")
brain_p1 <- LoadData("stxBrain", type = "posterior1")
brain <- merge(brain_a1, brain_p1)

# Calculating mitochondrial and hemoglobine genes percentage
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "mt.percent")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "hb.percent")

vlnplots <- VlnPlot(brain, features = c("nCount_Spatial", "nFeature_Spatial",
                                        "mt.percent", "hb.percent"),
                                        pt.size = 0.1, ncol = 2) + NoLegend()
vlnplots %>%
  ggsave(file = "vlnplots.svg", plot = ., width = 14, height = 14)

SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "mt.percent", "hb.percent")) %>%
  ggsave(file = "overview.svg", plot = ., width = 14, height = 28)

# Reading ncRNA features
ncrnas <- tolower(read.table("./Data/ncRNAs.txt", header=TRUE, sep="	", fill=TRUE)[,2])
ncrnas <- paste0(toupper(substr(ncrnas, 1, 1)), substr(ncrnas, 2, nchar(ncrnas)))

# Omitting empty
finding_relevant_feats <- function (brain_data, target_feats) {
  brain_counts <- brain_data[["Spatial"]]@counts
  feats <- intersect(rownames(brain_counts), target_feats)
  return(feats)
}
brain_ncrna_feats <- finding_relevant_feats(brain, ncrnas)

# Filtering data
brain <- brain[, brain$nFeature_Spatial > 500 & brain$mt.percent < 25 &
              brain$hb.percent < 20]

# Replotting data
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "mt.percent", "hb.percent")) %>%
  ggsave(file = "overview_filtered.svg", plot = ., width = 14, height = 28)


## Viewing genes with most counts
brain_counts <- brain@assays$Spatial@counts
brain_counts@x <- brain_counts@x/rep.int(colSums(brain_counts), diff(brain_counts@p))
most_expressed <- order(Matrix::rowSums(brain_counts), decreasing = T)[20:1]
svg("most_counts.svg")
boxplot(t(as.matrix(brain_counts[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

# Viewing ncRNAs with most counts
brain_ncrna <- subset(brain, features = brain_ncrna_feats)
ncrna_counts <- brain_ncrna@assays$Spatial@counts
ncrna_counts@x <- ncrna_counts@x/rep.int(colSums(ncrna_counts), diff(ncrna_counts@p))
most_expressed_ncrna <- order(Matrix::rowSums(ncrna_counts), decreasing = T)[20:1]
svg("ncrna_counts.svg")
ncrna_boxplot <- boxplot(t(as.matrix(ncrna_counts[most_expressed_ncrna, ])),
                         cex = 0.1, las = 1, xlab = "% total count per cell",
                         col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

# Filtering genes
brain <- brain[!grepl("Bc1", rownames(brain)), ]

brain <- brain[!grepl("^mt-", rownames(brain)), ]

norm_brain <- NormalizeData(brain)
# Finding variable features
norm_brain <- FindVariableFeatures(norm_brain, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(norm_brain), 10)
top10_ncrna <- intersect(VariableFeatures(norm_brain), brain_ncrna_feats)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(norm_brain)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
wrap_plots(plot1 + plot2) %>%
  ggsave(file = "variable_features.svg", plot = ., width = 14, height = 7)


SpatialFeaturePlot(brain, features = c("Hpca", "Ttr")) %>%
  ggsave(file = "spatial_features.svg", plot = ., width = 14, height = 7)



## Clustering integrated data
brain_list <- list(anterior1 = brain_a1, posterior1 = brain_p1)

# Running SCT on both datasets
brain_list <- lapply(brain_list, SCTransform, assay = "Spatial", method = "poisson")

# Setting maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # setting allowed size to 2K MiB


brain_features <- SelectIntegrationFeatures(brain_list, nfeatures = 3000, verbose = FALSE)
brain_list <- PrepSCTIntegration(object.list = brain_list, anchor.features = brain_features, verbose = FALSE)
brain_anchors <- FindIntegrationAnchors(object.list = brain_list, normalization.method = "SCT", verbose = FALSE, anchor.features = brain_features)
brain_integrated <- IntegrateData(anchorset = brain_anchors, normalization.method = "SCT", verbose = FALSE)

brain_integrated <- RunPCA(brain_integrated, verbose = FALSE)
brain_integrated <- FindNeighbors(brain_integrated, dims = 1:30)
brain_integrated <- FindClusters(brain_integrated, verbose = FALSE)
brain_integrated <- RunUMAP(brain_integrated, dims = 1:30)

DimPlot(brain_integrated, reduction = "umap", group.by = c("ident", "orig.ident")) %>%
  ggsave(file = "integrated_clusters.svg", plot = ., width = 14, height = 7)

SpatialDimPlot(brain_integrated) %>%
  ggsave(file = "integrated_clusters_spatial.svg", plot = ., width = 14, height = 7)

# Overlaying molecular data over tissue photo.
sorted_counts <- sort(rowSums(brain_integrated[["Spatial"]]@counts[brain_ncrna_feats,]), decreasing = TRUE)
brain_split_counts <- split(sorted_counts, ceiling(seq_along(sorted_counts)/4))
SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[1]])) %>%
    ggsave(file="ncRNA1.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[2]])) %>%
    ggsave(file="ncRNA2.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[3]])) %>%
    ggsave(file="ncRNA3.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[4]])) %>%
    ggsave(file="ncRNA4.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[5]])) %>%
    ggsave(file="ncRNA5.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[6]])) %>%
    ggsave(file="ncRNA6.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[7]])) %>%
    ggsave(file="ncRNA7.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[8]])) %>%
    ggsave(file="ncRNA8.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[9]])) %>%
    ggsave(file="ncRNA9.svg", plot=., width = 14, height = 28)

SpatialFeaturePlot(brain_integrated, features = names(brain_split_counts[[10]])) %>%
    ggsave(file="ncRNA10.svg", plot=., width = 14, height = 28)
