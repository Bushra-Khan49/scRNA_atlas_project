# Phase 3: Cell Type Annotation and Biomarker Validation

**Overview**: Having successfully mapped the dataset into distinct mathematical clusters, the final processing step is to assign true biological identities to these anonymous manifolds. After mapping the cellular lineages, I performed differential expression analysis to statistically validate the distinct transcriptomic signatures driving each cluster, ensuring rigorous scientific accuracy before final export.

---

## 1. Automated Biological Annotation

Historically, single-cell clusters were annotated manually by cross-referencing arbitrary gene lists—a process highly susceptible to human bias. To ensure robust, reproducible results, I utilized an automated reference-based annotation pipeline.

I employed the `SingleR` package to algorithmically compare the transcriptomic profile of my distinct clusters against the Human Primary Cell Atlas database. This effectively assigns highly accurate biological identities to the cells based on validated reference data rather than manual approximation.

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

# Add the predicted biological names back into the Seurat object
pbmc$cell_type <- pred$labels
```

---

## 2. Visualizing the Biological Cell Types

With the biological identities successfully mapped, I generated a definitive visualization of the cellular landscape. 

To ensure the visualization meets professional standards, I strictly maintained the legend format, keeping the text distinct from the graphical representation to prevent visual overlap on the smaller cellular populations.

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

## 3. Extracting Biomarker Signatures

Having identified the major cell populations, I next sought to statistically validate these annotations by uncovering their specific transcriptomic signatures. 

I executed a differential expression algorithm (`FindAllMarkers`) to test every single gene across all cell types. This identifies the specific biomarkers that are significantly upregulated in one cluster while remaining repressed in all others, generating a comprehensive statistical matrix of cellular identity.

```r
# Set the active identity to the new biological cell types
Idents(pbmc) <- "cell_type"

# Find all differentially expressed marker genes for every cell type
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the full table to a CSV file for rigorous downstream analysis
write.csv(markers, "../results/04_annotated_marker_genes.csv", row.names = FALSE)
```

---

## 4. Visualizing Transcriptomic Signatures (Heatmap)

To visually confirm the accuracy of the differential expression analysis, I extracted the top 5 most significantly upregulated genes for each cell type and plotted their expression profiles across the dataset.

Because clusters like T cells contain thousands of cells while others contain far fewer, a standard heatmap would visually misrepresent the smaller populations. To resolve this, I applied a strict downsampling parameter (`downsample = 30`), ensuring that every biological cell type is represented by an equal graphical width, standardizing the comparative visualization.

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

## 5. Final Export and Cloud Deployment Considerations

The final step is to save the fully annotated pipeline output. 

During the initial deployment of this project, I encountered a critical infrastructure issue: the standard `.rds` file was 274MB, which directly violated GitHub's 100MB strict file size limit, causing the cloud-hosted Shiny dashboard to crash.

To resolve this deployment failure, I utilized the `DietSeurat` function. This algorithmically strips the massive, intermediate assay matrices from the object while perfectly preserving the final cellular annotations and 2D UMAP coordinates. This successfully compressed the final object to ~10MB, allowing for seamless deployment.

```r
# Save the full annotated object for local analysis
saveRDS(pbmc, "../results/03_pbmc_annotated.rds")

# Shrink the object down to bypass GitHub's 100MB limit for cloud deployment
slim <- DietSeurat(pbmc, counts=TRUE, data=TRUE, scale.data=FALSE, dimreducs=c("umap"))
saveRDS(slim, "../results/03_pbmc_slim.rds")
```
