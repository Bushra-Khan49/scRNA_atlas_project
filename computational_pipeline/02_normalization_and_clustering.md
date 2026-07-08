This file takes the clean, filtered cells from the previous step and prepares them for grouping. I normalized the data so that the genes were actually comparable across different cells, found the genes that varied the most, and then grouped similar cells together into clusters.

---

# Phase 2: Normalization and Clustering

```r
library(Seurat)
library(ggplot2)
```

---

## Normalizing the Data
Different cells get sequenced at different depths just by chance, so comparing raw numbers directly doesn't work. I brought in the clean data from Phase 1 and normalized it. After that, I looked for the top 2,000 genes that varied the most between cells, because those are the ones that actually tell us what makes one cell different from another. I then scaled the data so that a few really loud genes wouldn't mess up the clustering.

```r
pbmc <- readRDS("../results/01_pbmc_filtered.rds")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
```

---

## Figuring out Principal Components (PCA)
Working with 2,000 genes is way too much for a computer to cluster easily. So, I ran PCA to squish the data down into the most important components.

**Troubleshooting Note**: Choosing how many PCs to use is always a bit tricky. If you use too few, you miss the biology. If you use too many, you just start grouping noise. I plotted an Elbow Plot to see where the variance dropped off. Looking at the graph, the line flattened out around 10, so I decided to stick with the first 10 PCs for my clustering. 

```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

ElbowPlot(pbmc, ndims = 30) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 3: PCA Elbow Plot", subtitle = "Variance explained by each Principal Component", x = "Principal Component", y = "Standard Deviation")
```

![Figure 3: PCA Variance (Elbow Plot)](../figures/02_elbow.png)


---

## Grouping the Cells (Clustering)
Once I had my 10 PCs, I used them to find which cells were closest to each other and grouped them into distinct clusters. 

To actually see these clusters, I ran UMAP, which takes all that complicated data and flattens it onto a simple 2D plot so I can visually check if the groups make sense.

```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

---

## What the Clusters Looked Like
**Figure 4: Mathematical Clustering**
When I plotted the UMAP, I got a nice visual showing distinct islands of cells. At this point, the program just grouped them by similarities but didn't know what kind of cells they actually were (it just called them Cluster 0, Cluster 1, etc.). I made sure to format the legend clearly so it was easy to read.

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

![Figure 4: Mathematical Clustering](../figures/02_umap_clusters.png)


---

## Saving the Clustered Data
I saved the object again. In the final step, I'll figure out what real cell types these clusters actually correspond to.

```r
saveRDS(pbmc, "../results/02_pbmc_clustered.rds")
```
