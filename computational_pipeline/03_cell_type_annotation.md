This file covers the final steps of the analysis. At this point, the cells were clustered into groups, but I didn't know what kind of cells they actually were. I used an automated tool to compare my clusters against a known database to figure out if they were T Cells, B Cells, etc., and then I found the specific genes that made each cell type unique.

---

# Phase 3: Cell Type Annotation and Finding Markers

```r
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(tidyverse)
```

---

## Figuring Out the Cell Types
Instead of manually guessing what cells were in each cluster by looking at genes one by one (which takes forever and is easy to mess up), I used a package called `SingleR`. It takes my cells and checks them against the Human Primary Cell Atlas database to automatically assign the correct names.

```r
pbmc <- readRDS("../results/02_pbmc_clustered.rds")

ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(pbmc)

pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
pbmc$cell_type <- pred$labels
```

---

## What the Final Labeled Map Looked Like
**Figure 5: Biological Cell Type Map**
After running the annotation, I mapped the real cell names back onto my UMAP plot. 

**Troubleshooting Note**: When I first made this plot, Seurat tried to put the text labels directly on top of the dots. It looked incredibly messy and overlapping, and I couldn't read the smaller clusters. I fixed it by turning the labels off on the plot itself and instead moving them to a nice, bold legend on the side so the actual data points are clearly visible.

```r
DimPlot(pbmc, group.by = "cell_type", label = FALSE, pt.size = 0.5) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(title = "Figure 5: Biological Cell Type Map", subtitle = "UMAP projection mapped via SingleR algorithmic annotation", x = "UMAP Dimension 1", y = "UMAP Dimension 2")
```

![Figure 5: Biological Cell Type Map](../figures/03_umap_celltype.png)


---

## Finding the Marker Genes
Now that I knew what the cells were, I wanted to find the exact genes that were driving those identities. I ran a test across all the cell types to find the genes that were highly turned on in one specific cell type but turned off everywhere else. I saved these genes into a CSV file.

```r
Idents(pbmc) <- "cell_type"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "../results/04_annotated_marker_genes.csv", row.names = FALSE)
```

---

## Visualizing the Marker Genes
**Figure 6: Biomarker Expression Heatmap**
To double check that the markers were right, I grabbed the top 5 genes for each cell type and plotted them on a heatmap. 

**Troubleshooting Note**: I originally tried plotting the heatmap *before* the cells were annotated, which just gave me a legend with meaningless numbers like 0, 1, and 2. Also, because some clusters had tons of cells (like T Cells) and others had very few (like Platelets), the columns for the small clusters got completely squished and were impossible to see. I fixed this by running the heatmap *after* annotation so the real names show up, and I forced the plot to only sample exactly 30 cells from each group so the columns all have the exact same width.

```r
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

pbmc_downsampled <- subset(pbmc, downsample = 30)

DoHeatmap(pbmc_downsampled, features = top5$gene, size = 4, angle = 90) +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 6: Biomarker Expression Heatmap", subtitle = "Top 5 differentially expressed genes per biological cell type")
```

![Figure 6: Biomarker Expression Heatmap](../figures/03_marker_heatmap.png)


---

## Final Save and Dashboard Fix
I saved the final fully-annotated data object.

**Troubleshooting Note**: This final `.rds` file ended up being 274MB. When I tried to push this to GitHub to host my Shiny dashboard, it completely failed because GitHub has a strict 100MB file limit. Because the file didn't upload, my dashboard couldn't find the data and crashed. I fixed this by using the `DietSeurat` function to strip out all the heavy background matrices I didn't need for the dashboard, shrinking the file down to a tiny 10MB while keeping the cell names and UMAP coordinates intact.

```r
saveRDS(pbmc, "../results/03_pbmc_annotated.rds")

slim <- DietSeurat(pbmc, counts=TRUE, data=TRUE, scale.data=FALSE, dimreducs=c("umap"))
saveRDS(slim, "../results/03_pbmc_slim.rds")
```
