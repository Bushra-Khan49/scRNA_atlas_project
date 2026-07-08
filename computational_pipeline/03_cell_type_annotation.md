# Phase 3: Cell Type Annotation and Finding Markers

**Overview**: In this final processing file, I am taking the anonymous, mathematically generated clusters from Phase 2 and figuring out what actual cells they are (T Cells, B Cells, Monocytes, etc.). After naming the clusters, I run an algorithm to pull out the exact biomarker genes that make each of those cell types unique.

---

## 1. Automated Biological Annotation

In the old days, scientists used to guess what the clusters were by manually staring at gene lists. This is a terrible idea because it is incredibly prone to human bias. 

Instead, I use the `SingleR` package. Make sure you do this: `SingleR` takes your data and automatically compares it against a massively curated database (the Human Primary Cell Atlas). It checks the transcriptomic profile of your cells against the known database and statistically predicts the true biological names of your clusters. 

```r
# Load required libraries
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(tidyverse)

# Load the clustered data from Phase 2
pbmc <- readRDS("../results/02_pbmc_clustered.rds")

# Pull down the Human Primary Cell Atlas reference database
ref <- celldex::HumanPrimaryCellAtlasData()
# Convert Seurat object to SingleCellExperiment format required by SingleR
sce <- as.SingleCellExperiment(pbmc)

# Run the automated annotation algorithm
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# Add the predicted biological names back into our Seurat object
pbmc$cell_type <- pred$labels
```

---

## 2. Visualizing the Biological Cell Types

Now that the cells have real names, I re-plot the UMAP. 

Make sure you pay attention to the legend formatting here. If you just slap the text labels on the plot, it becomes a jumbled mess where you can't even read the names of the small cell populations. I moved the labels to a bold side legend so that the actual visual representation of the cells remains completely clean.

```r
# Plot the UMAP colored by the new biological cell types
DimPlot(pbmc, group.by = "cell_type", label = FALSE, pt.size = 0.5) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(title = "Figure 5: Biological Cell Type Map", subtitle = "UMAP projection mapped via SingleR algorithmic annotation", x = "UMAP Dimension 1", y = "UMAP Dimension 2")
```

![Figure 5: Biological Cell Type Map](../figures/03_umap_celltype.png)

---

## 3. Finding the Biomarker Genes

Now that I know what the cells are, I need to know *why* they are what they are. I run `FindAllMarkers`, which runs a statistical test across every single cell type. It looks for genes that are highly expressed in one group but completely turned off everywhere else. I then save this massive table of genes as a CSV file for downstream use.

```r
# Set the active identity to our new biological cell types
Idents(pbmc) <- "cell_type"

# Find all differentially expressed marker genes for every cell type
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the full table to a CSV file
write.csv(markers, "../results/04_annotated_marker_genes.csv", row.names = FALSE)
```

---

## 4. Visualizing the Biomarkers (Heatmap)

To visually prove that my cell types are correct, I grab the top 5 most highly expressed genes for each group and plot them on a heatmap. 

You must be careful here: because some cell types (like T Cells) have thousands of cells and others have very few, a normal heatmap will completely crush the small clusters out of existence. I fixed this by using the `downsample = 30` parameter, which forces the plot to sample exactly 30 cells from each group. This guarantees that every single cell type has the exact same visual width on the graph.

```r
# Extract the top 5 marker genes per cluster
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

# Downsample the cells so the heatmap columns are perfectly balanced
pbmc_downsampled <- subset(pbmc, downsample = 30)

# Plot the heatmap
DoHeatmap(pbmc_downsampled, features = top5$gene, size = 4, angle = 90) +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 6: Biomarker Expression Heatmap", subtitle = "Top 5 differentially expressed genes per biological cell type")
```

![Figure 6: Biomarker Expression Heatmap](../figures/03_marker_heatmap.png)

---

## 5. Final Export and Cloud Dashboard Fix

Finally, I save the fully annotated object. 

I actually had to do some serious troubleshooting here. My main `.rds` file was 274MB. When I tried to push this to GitHub to host my interactive Shiny dashboard, GitHub blocked it because of their 100MB file limit. Because it wouldn't upload, my dashboard completely crashed. 

To fix this, you should use the `DietSeurat` function. It strips out all the massive background math matrices that we no longer need, shrinking the file down to just 10MB while keeping the cell names and UMAP coordinates perfectly intact for the dashboard.

```r
# Save the full annotated object for local use
saveRDS(pbmc, "../results/03_pbmc_annotated.rds")

# Shrink the object down to bypass GitHub's 100MB limit for cloud deployment
slim <- DietSeurat(pbmc, counts=TRUE, data=TRUE, scale.data=FALSE, dimreducs=c("umap"))
saveRDS(slim, "../results/03_pbmc_slim.rds")
```
