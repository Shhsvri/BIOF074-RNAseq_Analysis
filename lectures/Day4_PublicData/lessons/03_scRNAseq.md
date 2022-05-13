---
title: single-cell RNA sequencing analysis using Seurat
author: Shahin Shahsavari
date: 05/09/2022
duration: 120 minutes 
---


### Analysis of scRNA-seq in R using Seurat

For RNAseq data normalization and analysis, we will be using the
[Seurat](https://satijalab.org/seurat/) library.

### Load the required libraries

```R
# load in the libraries
library(tidyverse)
library(Seurat)
library(patchwork) # an extension of ggplot2
```


### Load the PBMC dataset

Similar to DESeq2, Seurat creates a seurat-specific object for storage of raw counts, metadata,
transformation,, and further downstream analysis

```R
setwd("~/Day4")
pbmc_data <- Read10X(data.dir = "scRNAseq_data/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

To get the names of all the slots in a Seurat object, use `slotNames()`

```R
# Inspect the Seurat object using slotNames()
slotNames(pbmc)
# access the metadata and assays using the `@` operator
pbmc@meta.data %>% head()
pbmc@assays
slotNames(pbmc@assays$RNA)
```

Let's examine the data

```R
# we subset a few genes of interest in the first 5 cells
pbmc_data[c("CD3D", "TCL1A", "MS4A1"), 1:5]
```
The `.` values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are
0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed
savings for Drop-seq/inDrop/10x data.

### Standard Pre-processing and Quality Control

Seurat allows you to explore QC metrics and filter cells based on a variety of criteria. The common ones are:

* The number of unique genes detected in a cell `nFeature_RNA`
	* low-quality cells or empty droplets tend to have very few genes
	* doublets or multiplets often have an abnormally high gene count
* The total number of molecules within a cell `nCount_RNA`
	* number of features (genes) correlates strongly with unique genes
* The percentage of reads that map to the mitochondrial genome
	* high mitochondrial RNA in a sample signify cell death (>5-10%)

Let's add the mitochondrial RNA data to our Seurat object.
Then we could visualize the QC metrics using violin plots

```R
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

Alternatively, you could visualize the QC metrics using `FeatureScatter`:

```R
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

We extract cells that have feature counts between 200 and 2,500

```R
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

### Data Normalization

After QC, we normalize our raw count matrix using `NormalizeData`. This uses a global-scaling
normalization method called "LogNormalize".

```R
pbmc <- NormalizeData(pbmc)
```

### Identification of highly variable features (genes)

Below, we extract the genes that have the highest normalized counts variable between all cells. In order to speed
up this process, this returns 2,000 genes per dataset. These will be used further downstream during our analysis.

```R
pbmc <- FindVariableFeatures(pbmc,
				   selection.method = "vst",
				   nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

### Center scaling the data

Similar to bulkRNAseq, we need to center our values around 0 and the variance to 1 for all genes. This is required
for heatmaps and clustering.

```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc)
```

### Perform linear dimensionality reduction

Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input,
but can be defined using features argument if you wish to choose a different subset.

```R
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Visualize the standard deviation along the top 20 PCs
ElbowPlot(pbmc)
```

### Clustering

Clustering is one of the core steps of scRNA-seq analysis. Using our reduced dimensions, want to find subpopulations within
our set of data that correspond to distinct cell types, which we do by clustering. Seurat employs the Louvain clustering
algorithm for clustering.

```R
# Use the first 10 PCs for clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# resolution of 0.5 is optimized for ~ 3000 single-cell datasets
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```

### Run non-linear dimensionality reduction

Other dimensionality reduction techniques such as tSNE (t-stochastic neighbor embedding) or UMAP (Uniform Manifold Approximation and
Projection) are offered through Seurat. UMAP has been the most common method for scRNAseq data visualization.

```R
# Add the UMAP values to our pbmc object
pbmc <- RunUMAP(pbmc, dims = 1:10)# only first 10
# Plot UMAP
DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident')
# Now add the clusters
DimPlot(pbmc, reduction = "umap", group.by ='seurat_clusters' )
# Identifying markers for each cluster. these are the genes that are representative of only 1 cluster.
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers)

# Generate violin plots, ridge plots, and feature plots of these clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), group.by ='seurat_clusters' )
RidgePlot(pbmc, features = c("NKG7", "PF4"),group.by ='seurat_clusters' )
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
    "CD8A"))

top10 <- pbmc.markers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC) # select top 10 markers for each gene
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```


