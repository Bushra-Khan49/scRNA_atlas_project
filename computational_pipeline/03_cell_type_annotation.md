# Phase 3: Biological Annotation and Biomarker Discovery

**Overview**: In Phase 2, the cells were grouped mathematically into anonymous clusters based on transcriptomic similarity. In this final phase, my goal was to assign real biological identities (e.g., T Cells, B Cells, Monocytes) to these groups, and to extract the specific biomarker genes driving those identities.

```r
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(tidyverse)
```

## Algorithmic Biological Annotation
Historically, annotating single-cell data involved manually running differential expression, staring at marker genes, and subjectively assigning cell types. This is highly prone to human error and bias.

To ensure strict scientific reproducibility, I used the `SingleR` algorithmic pipeline. `SingleR` compares the transcriptomic profile of my cells against the deeply curated Human Primary Cell Atlas (HPCA) database, automatically assigning statistically significant biological labels.

```r
pbmc <- readRDS("../results/02_pbmc_clustered.rds")

# Pull the curated reference database
ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(pbmc)

# Run the algorithmic prediction
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
pbmc$cell_type <- pred$labels
```

## Results & Interpretation: Biologically Annotated UMAP
**Figure 5: Biological Cell Type Map**
I mapped the algorithm's predictions back onto my original Seurat object. 

**Visual Fix**: When I initially plotted this, having the text labels directly on top of the cell clusters made the diagram look messy and unreadable. I completely removed the on-plot text and built a bold, highly readable legend on the right side. This preserves the purity of the data points while remaining perfectly decipherable. We can clearly observe the massive CD4+ and CD8+ T Cell populations, along with distinct islands of B Cells and Monocytes.

```r
DimPlot(pbmc, group.by = "cell_type", label = FALSE, pt.size = 0.5) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(title = "Figure 5: Biological Cell Type Map", subtitle = "UMAP projection mapped via SingleR algorithmic annotation", x = "UMAP Dimension 1", y = "UMAP Dimension 2")
```

## Biomarker Discovery (Differential Expression)
Now that the cells have real biological names, I ran a Wilcoxon Rank Sum test across the annotated groups to find the "biomarkers"—the specific genes that are highly expressed in one cell type but turned off in all the others.

```r
Idents(pbmc) <- "cell_type"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "../results/04_annotated_marker_genes.csv", row.names = FALSE)
```

## Results & Interpretation: Biomarker Expression Heatmap
**Figure 6: Biomarker Expression Heatmap**
To prove the differential expression was successful, I extracted the top 5 most defining genes for each biological cell type and plotted them.

**Visual Fix**: Previously, I plotted the heatmap *before* annotation, resulting in a meaningless legend (0, 1, 2). Furthermore, rare cell types (like Platelets) had their columns physically crushed out of existence by massive populations (like T Cells). 
By moving the heatmap to *after* the biological annotation, and applying a strict `downsample = 30` parameter, I guaranteed that every single cell type receives the exact same column width. The biological names are bolded and perfectly readable. The bright yellow blocks definitively prove the specificity of our top 5 biomarkers per cell type.

```r
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

# Downsample to exactly 30 cells per biological identity to ensure perfect visual balance
pbmc_downsampled <- subset(pbmc, downsample = 30)

DoHeatmap(pbmc_downsampled, features = top5$gene, size = 4, angle = 90) +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 6: Biomarker Expression Heatmap", subtitle = "Top 5 differentially expressed genes per biological cell type")
```

## Final Data Checkpoint
I saved the fully processed object. 
**Deployment Optimization**: I generated both the full 274MB version and a highly compressed "slim" 10MB version (stripping out raw matrices and preserving only UMAP coordinates/identities). This was necessary because the 274MB file breached GitHub's 100MB limit, causing an out-of-memory crash when deploying the final Shiny Dashboard to Posit Connect Cloud.

```r
saveRDS(pbmc, "../results/03_pbmc_annotated.rds")

# Generate slim version for Cloud deployment
slim <- DietSeurat(pbmc, counts=TRUE, data=TRUE, scale.data=FALSE, dimreducs=c("umap"))
saveRDS(slim, "../results/03_pbmc_slim.rds")
```
