# Phase 1: Quality Control

**Overview**: Single-cell RNA sequencing data is inherently noisy. Before conducting any downstream clustering, we must rigorously filter out droplets containing no real cells (empty droplets), droplets containing multiple cells (doublets), and dying cells whose membranes have ruptured.

## Environment Setup
Setting up the environment for scRNA-seq analysis in R is notoriously tricky due to strict Bioconductor versioning. 
**Troubleshooting Note**: When installing dependencies, `SingleCellExperiment` continuously caused crashes during the annotation phase because the default CRAN version was out of date. I resolved this by explicitly forcing the Bioconductor version via `BiocManager::install("SingleCellExperiment")`. Always ensure your `renv` is synced!

```r
library(Seurat)
library(tidyverse)
library(ggplot2)
```

## Loading the Raw Data
I loaded the raw 10x Genomics PBMC dataset and initialized a Seurat object. Right at initialization, I applied two extremely lenient baseline filters to immediately discard massive amounts of noise:
- `min.cells = 3`: Removes any gene that isn't expressed in at least 3 cells across the entire dataset.
- `min.features = 200`: Removes any cell with fewer than 200 distinct genes detected (almost certainly empty fluid droplets).

```r
mat_path <- "../data/pbmc3k/filtered_gene_bc_matrices/hg19"
counts   <- Read10X(data.dir = mat_path)

pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)
```

## Mitochondrial Percentage (Cell Health)
When a cell is stressed or dying, its membrane ruptures and cytoplasmic RNA leaks out, but the mitochondria remain trapped inside. Thus, an unusually high ratio of mitochondrial RNA to total RNA is the standard biomarker for a dead cell.
I calculated this percentage by searching for all genes starting with `MT-`.

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

## Results & Interpretation: Pre-Filter Quality Assessment
Before applying strict cutoffs, I plotted the raw distributions.

**Figure 1: Quality Control Metrics**
This violin plot shows the spread of features (unique genes), counts (total RNA molecules), and the mitochondrial percentage per cell. We can clearly see a dense cluster of healthy cells, but a long "tail" of dying cells with >10% MT.

```r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, cols = "#3A86FF") +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 1: Quality Control Metrics", subtitle = "Distribution of RNA counts and Mitochondrial percentage", y = "Count / Percentage")
```

**Figure 2: Mitochondrial vs Total RNA**
This scatter plot maps the total RNA counts against the MT percentage. The cells in the upper-left quadrant (low RNA, high MT) are definitively dead and must be removed.

```r
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = "#FF006E") +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 2: Mitochondrial vs Total RNA", subtitle = "Scatter plot identifying stressed/dying cells (high MT%)", x = "Total RNA Molecules Detected", y = "Mitochondrial Percentage (%)")
```

## Applying Biological Cutoffs
Based on the visual evidence above, I applied the following strict biological cutoffs:
1. `nFeature_RNA > 200`: Confirming we only keep true cells.
2. `nFeature_RNA < 2500`: Filtering out doublets (two cells trapped in one droplet will have an impossibly high gene count).
3. `percent.mt < 5`: Dropping the dead cells identified in Figure 2.

```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

## Saving Checkpoint
I saved the filtered matrix. In large-scale single-cell pipelines, you must save checkpoints after QC to avoid re-running intensive filtering algorithms every time you tweak a downstream parameter.

```r
saveRDS(pbmc, "../results/01_pbmc_filtered.rds")
```
