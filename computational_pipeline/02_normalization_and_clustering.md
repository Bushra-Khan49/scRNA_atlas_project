# Phase 2: Normalization and Clustering

**Overview**: In this file, I am taking the clean, healthy cells that survived our quality control phase, and figuring out how to group them together. To do this, you first have to normalize the data so you can actually compare the cells to each other fairly. Then, you pull out the specific genes that make the cells different from one another, crush the data down to its most important dimensions using PCA, and finally cluster them into distinct biological islands.

---

## 1. Normalizing the Data

Whenever you run single-cell sequencing, different cells end up getting sequenced at slightly different depths purely by chance. Some might have 1,000 RNA molecules detected, and others might have 1,500, even if they are identical cells. Make sure you keep in mind that if you don't normalize this, your computer will think those cells are different just because of a technical error. 

Here, I load the clean data and apply a standard Log Normalization.

```r
# Load Seurat and ggplot2
library(Seurat)
library(ggplot2)

# Load the filtered data object we saved in Phase 1
pbmc <- readRDS("../results/01_pbmc_filtered.rds")

# Normalize the data so gene expression is comparable across all cells
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

---

## 2. Finding the Variable Genes

Out of the 20,000+ genes in the human genome, most of them do the exact same basic housekeeping tasks in every single cell. We don't care about those. We only want the genes that are highly turned on in some cells but turned off in others. 

I tell Seurat to find the top 2,000 most highly variable genes. Then, I scale the data. You must do this, because otherwise, a gene that is naturally expressed at very high levels will completely overpower a gene that is expressed at low levels but is actually very biologically important.

```r
# Identify the top 2000 genes with the highest variance between cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data so the mean expression is 0 and variance is 1
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
```

---

## 3. Dimensionality Reduction (PCA)

Trying to cluster cells across 2,000 different gene dimensions is nearly impossible for a computer. It's called the "curse of dimensionality". To fix this, I run Principal Component Analysis (PCA), which squishes all that variance down into a few super-components.

How do you know how many Principal Components to use? You run an Elbow Plot. 

```r
# Run PCA on the scaled data using our 2000 variable features
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

Looking at the graph above, the variance drops sharply and then flattens out around 10. That's why I chose to use the first 10 PCs for the final clustering. 

---

## 4. Grouping the Cells (Clustering)

Now that I have my 10 dimensions, I build a K-Nearest Neighbor graph, which essentially draws lines between the cells that are most similar to each other. Then I run a clustering algorithm to group them into actual communities.

Finally, I run UMAP. UMAP takes that 10-dimensional graph and crushes it down to a 2D map so human eyes can actually look at it and verify the clusters.

```r
# Find the closest neighbors in 10-dimensional PCA space
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Group the neighbors into distinct clusters (resolution 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP to project the clusters onto a 2D plot
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

---

## 5. Visualizing the Clusters

Here, I generate the UMAP plot. At this stage, the computer doesn't know what a "T Cell" or a "B Cell" is; it just knows that certain cells look mathematically similar. That's why the legend just says "Cluster 0", "Cluster 1", etc. 

Notice that I turned the text labels *off* the plot and put them in a legend instead. You should always do this if you have a lot of clusters, otherwise the text overlaps the dots and looks extremely messy.

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

## 6. Saving the Clustered Data

I saved this clustered object as a new file. In the final phase, we are going to figure out what those clusters actually represent biologically!

```r
# Save the clustered Seurat object
saveRDS(pbmc, "../results/02_pbmc_clustered.rds")
```
