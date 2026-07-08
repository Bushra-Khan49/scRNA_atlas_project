# Phase 3: Cell Type Annotation and Biomarker Validation

**Project Context**: In Phase 2, we successfully separated the blood sample into distinct mathematical clusters based on transcriptomic variance. However, "Cluster 0" and "Cluster 1" mean nothing in a biological context. The final objective of this pipeline is to construct a biologically accurate map of the immune system. Therefore, this phase takes those anonymous mathematical clusters, algorithmically identifies their true biological lineages (e.g., T Cells, B Cells), and statistically validates those identities by extracting their defining biomarker genes. 

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

```r
# Generate the UMAP colored by the new biological cell types
umap_celltype <- DimPlot(pbmc, group.by = "cell_type", label = FALSE, pt.size = 0.5) +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(title = "Figure 5: Biological Cell Type Map", subtitle = "UMAP projection mapped via SingleR algorithmic annotation", x = "UMAP Dimension 1", y = "UMAP Dimension 2")

ggsave("../figures/03_umap_celltype.png", plot = umap_celltype, width = 9, height = 6, bg = "white")
```

![Figure 5: Biological Cell Type Map](../figures/03_umap_celltype.png)

### Image Interpretation: Identifying the Lineages
Looking at **Figure 5**, you can see that the algorithm successfully assigned a biological identity to every single island we found in Phase 2. The massive teal cluster dominating the left side of the map represents T Cells. On the right, we see perfectly isolated pockets of B Cells (orange), Monocytes (purple), and NK Cells (blue). Because these colors perfectly map to the physical islands we saw previously, it visually confirms that our mathematically distinct clusters actually represent entirely different cell lineages in the human immune system.

---

## 3. Extracting Biomarker Signatures

Having identified the major cell populations visually, I next sought to statistically validate these annotations by uncovering their specific transcriptomic signatures. 

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

Because clusters like T cells contain thousands of cells while others contain far fewer, a standard heatmap would visually crush the smaller populations. I fixed this by applying a strict downsampling parameter (`downsample = 30`), ensuring that every biological cell type has the exact same graphical width.

```r
# Extract the top 5 marker genes per cluster
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

# Downsample the cells so the heatmap columns are perfectly balanced
pbmc_downsampled <- subset(pbmc, downsample = 30)

# Generate the heatmap
marker_heatmap <- DoHeatmap(pbmc_downsampled, features = top5$gene, size = 4, angle = 90) +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 6: Biomarker Expression Heatmap", subtitle = "Top 5 differentially expressed genes per biological cell type")

ggsave("../figures/03_marker_heatmap.png", plot = marker_heatmap, width = 10, height = 7, bg = "white")
```

![Figure 6: Biomarker Expression Heatmap](../figures/03_marker_heatmap.png)

### Image Interpretation: Validating the Cell Types
Looking at **Figure 6**, the bright yellow color represents extremely high gene expression, while purple represents no expression. Notice how the yellow forms distinct, solid blocks that perfectly align with the columns (the cell types). If the annotation algorithm was wrong, the yellow would be scattered randomly everywhere. 

For instance, if you look at the B_cell column on the far right, you see a massive block of bright yellow specifically for genes like CD79A. Because CD79A is a universal, biological requirement for B cell survival, seeing it physically concentrated *only* in the B cell column provides visual, definitive proof that our cell types are correctly annotated.

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
