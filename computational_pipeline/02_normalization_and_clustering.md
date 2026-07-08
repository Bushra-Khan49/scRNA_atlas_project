# Phase 2: Normalization and Clustering

**Overview**: Having isolated a highly viable cellular population during the quality control phase, I now proceed to group these cells based on their transcriptomic similarities. This requires normalizing the raw expression data to account for technical sequencing variances, identifying the genes that drive biological diversity, and performing dimensionality reduction to map the cells into distinct mathematical clusters.

---

## 1. Normalization and Variance Stabilization

In single-cell sequencing, individual cells are captured and sequenced at varying depths due to technical artifacts. Comparing raw expression counts directly would falsely identify deeply sequenced cells as biologically distinct from shallowly sequenced cells. 

To correct for this, I applied Log Normalization to the filtered dataset. Following normalization, I isolated the top 2,000 highly variable genes. Because most of the ~20,000 genes in the genome are ubiquitously expressed housekeeping genes, focusing exclusively on these 2,000 variable features allows the algorithm to cluster cells based purely on true biological variance. Finally, the data was scaled so that highly expressed genes do not artificially dominate the clustering algorithm.

```r
# Load Seurat and ggplot2
library(Seurat)
library(ggplot2)

# Load the filtered data object saved in Phase 1
pbmc <- readRDS("../results/01_pbmc_filtered.rds")

# Normalize the data so gene expression is comparable across all cells
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the top 2000 genes with the highest variance between cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data so the mean expression is 0 and variance is 1
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
```

---

## 2. Dimensionality Reduction (PCA)

Attempting to cluster cells across 2,000 dimensions simultaneously is computationally inefficient and highly susceptible to noise. I utilized Principal Component Analysis (PCA) to compress this variance into orthogonal principal components.

To determine the optimal number of dimensions for downstream clustering, I generated an Elbow Plot to visualize the variance captured by each principal component.

```r
# Run PCA on the scaled data using the 2000 variable features
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Plot the variance explained by each Principal Component
ElbowPlot(pbmc, ndims = 30) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", size = 1) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 3: PCA Variance (Elbow Plot)", subtitle = "Dashed red line indicates dimensional cutoff at PC 10", x = "Principal Component Index (Dimensions)", y = "Standard Deviation (Variance)")
```

![Figure 3: PCA Variance (Elbow Plot)](../figures/02_elbow.png)

As demonstrated above, the standard deviation drops sharply and begins to plateau around PC 10. Therefore, I selected the first 10 dimensions to capture the majority of the true biological signal while excluding downstream statistical noise.

---

## 3. Unsupervised Clustering

Utilizing these 10 dimensions, I constructed a K-Nearest Neighbor (KNN) graph to identify cells with similar transcriptomic profiles. A Louvain clustering algorithm was then applied to partition this graph into distinct cellular communities.

To visualize these multidimensional clusters, I projected the data into a 2D space using Uniform Manifold Approximation and Projection (UMAP).

```r
# Find the closest neighbors in 10-dimensional PCA space
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Group the neighbors into distinct clusters (resolution 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP to project the clusters onto a 2D plot
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

---

## 4. Visualizing the Mathematical Clusters

At this stage, the clustering is purely mathematical. The algorithm has grouped the cells by transcriptomic similarity, assigning them anonymous labels (e.g., "Cluster 0"). 

To ensure the visualization remains clean and readable, I configured the plot to display a distinct legend rather than overlaying text directly onto the data points, which can obscure small cellular populations.

```r
# Update the cluster names so they look cleaner in the legend
pbmc$seurat_clusters <- paste("Cluster", pbmc$seurat_clusters)
Idents(pbmc) <- "seurat_clusters"

# Plot the UMAP
DimPlot(pbmc, reduction = "umap", label = FALSE, pt.size = 0.5) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 10)) +
  labs(title = "Figure 4: Mathematical Clustering", subtitle = "UMAP projection of KNN louvain clusters", x = "UMAP Dimension 1", y = "UMAP Dimension 2")
```

![Figure 4: Mathematical Clustering](../figures/02_umap_clusters.png)

---

## 5. Saving the Clustered Dataset

I saved the clustered Seurat object. In the final phase of the pipeline, I will transition from mathematical clustering to biological annotation, identifying exactly what cell types these distinct islands represent.

```r
# Save the clustered Seurat object
saveRDS(pbmc, "../results/02_pbmc_clustered.rds")
```
