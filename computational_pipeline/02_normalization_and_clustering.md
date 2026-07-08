# Phase 2: Mathematical Normalization and Clustering

**Overview**: After establishing a clean dataset of viable cells, the next stage of the pipeline transforms the raw RNA sequencing counts mathematically. The goal is to normalize technical biases, identify the most critical genes, compress the dimensionality of the data, and physically cluster the cells based on their transcriptomic similarities.

```r
library(Seurat)
library(ggplot2)
```

## Normalization and Variance Scaling
Raw sequencing data is heavily biased by sequencing depth—some cells are randomly sequenced much deeper than others. 
I normalized the data using `LogNormalize` (scaling to 10,000 to match standard TPM rates).

```r
pbmc <- readRDS("../results/01_pbmc_filtered.rds")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

Out of the ~20,000 genes in the human genome, most are just biological background noise (e.g., housekeeping genes). I ran a Variance Stabilizing Transformation (VST) to extract the top 2,000 most highly variable genes. These are the critical genes driving the biological differences between cell types. I then scaled the data (mean = 0, variance = 1) so that highly expressed genes wouldn't artificially dominate the algorithm.

```r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
```

## Dimensionality Reduction (PCA)
A 2000-dimension matrix (one dimension per variable gene) is computationally impossible to cluster effectively due to the curse of dimensionality. I ran Principal Component Analysis (PCA) to compress the data down to linear components of maximum variance.

```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
```

## Results & Interpretation: PCA Variance
**Figure 3: PCA Elbow Plot**
To determine how many Principal Components to keep for clustering, I generated an Elbow Plot. The y-axis shows the standard deviation (variance explained) by each PC. We observe a sharp "elbow" drop-off around PC 10. Therefore, I determined that the vast majority of true biological signal is captured in the first 10 PCs, and anything beyond that is largely noise.

```r
ElbowPlot(pbmc, ndims = 30) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 3: PCA Elbow Plot", subtitle = "Variance explained by each Principal Component", x = "Principal Component", y = "Standard Deviation")
```

## Manifold Approximation and Clustering
Using those top 10 PCs, I mapped the cells in multidimensional space by building a K-Nearest Neighbor (KNN) graph. I then applied the Louvain community detection algorithm (`resolution = 0.5`) to group these neighbors into distinct, anonymous mathematical clusters.

To visualize these clusters, I applied Uniform Manifold Approximation and Projection (UMAP), a non-linear dimensionality reduction technique that crushes the 10D space down to a 2D plane.

```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

## Results & Interpretation: Mathematical Clustering
**Figure 4: Mathematical Clustering**
This UMAP projection reveals distinct "islands" of cells. At this stage, the algorithm has successfully grouped the cells by their transcriptomic profiles, but it does not know what they are. I updated the metadata so the legend correctly reads "Cluster X" rather than raw, undefined numbers, and utilized `ggrepel` to ensure the bold text labels point directly to the centroids without obscuring the data.

```r
pbmc$seurat_clusters <- paste("Cluster", pbmc$seurat_clusters)
Idents(pbmc) <- "seurat_clusters"

DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, pt.size = 0.5) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold")) +
  labs(title = "Figure 4: Mathematical Clustering", subtitle = "UMAP projection of K-Nearest Neighbor louvain clusters", x = "UMAP Dimension 1", y = "UMAP Dimension 2")
```

## Saving Checkpoint
With the cells successfully clustered, I saved the object to pass to the final biological annotation pipeline.

```r
saveRDS(pbmc, "../results/02_pbmc_clustered.rds")
```
