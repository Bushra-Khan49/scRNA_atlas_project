This file covers the very first part of my analysis. I started by getting the raw single-cell data, loading it into R, and doing some basic quality checks to filter out bad or dying cells before doing any real analysis.

---

# Phase 1: Quality Control

Before doing anything, I set up my environment. I loaded Seurat (which is basically the standard for this kind of data) and the tidyverse package to help with making plots and handling data.

```r
library(Seurat)
library(tidyverse)
library(ggplot2)
```

---

## Loading the Raw Data
I downloaded the raw dataset from 10x Genomics. The files came as raw matrices, which I then loaded into R. Since single-cell data has a lot of empty droplets (which happen during the sequencing process), I told Seurat to immediately drop anything that didn't look like a real cell just to save memory. 

```r
mat_path <- "../data/pbmc3k/filtered_gene_bc_matrices/hg19"
counts   <- Read10X(data.dir = mat_path)

pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)
```

---

## Checking for Dead Cells (Mitochondrial RNA)
One of the biggest issues with this kind of data is that dying cells leak their regular RNA but keep their mitochondrial RNA. So, if a cell has a high percentage of mitochondrial genes, it's usually dead or dying. I calculated this percentage by looking for genes that start with "MT-".

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

---

## What the Data Looked Like Before Filtering
I wanted to check what the data actually looked like before I started dropping cells. I made some plots to see the distribution of the genes and the mitochondrial percentage.

**Figure 1: Quality Control Metrics**
This plot showed me that most cells were healthy and grouped together, but there was a clear group of dying cells with really high mitochondrial percentages that I needed to get rid of.

```r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, cols = "#3A86FF") +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 1: Quality Control Metrics", subtitle = "Distribution of RNA counts and Mitochondrial percentage", y = "Count / Percentage")
```

**Figure 2: Mitochondrial vs Total RNA**
This scatter plot confirmed the same thing. The dots in the upper left corner had hardly any real RNA but tons of mitochondrial RNA, meaning they were definitely dead.

```r
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = "#FF006E") +
  theme_classic() +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 2: Mitochondrial vs Total RNA", subtitle = "Scatter plot identifying stressed/dying cells (high MT%)", x = "Total RNA Molecules Detected", y = "Mitochondrial Percentage (%)")
```

---

## Filtering Out the Bad Data
**Troubleshooting Note**: Figuring out the exact cutoff numbers was a bit tricky. If I was too strict, I'd lose too much good data. If I was too loose, the downstream analysis would be noisy. Based on the plots, I decided to filter out anything with more than 5% mitochondrial RNA, and dropped cells with too many or too few genes to avoid doublets (two cells stuck together) and empty droplets.

```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

---

## Saving the Clean Data
Since processing the data takes time, I saved this clean version so I wouldn't have to keep re-running the filtering steps every time I wanted to test something in the next phase.

```r
saveRDS(pbmc, "../results/01_pbmc_filtered.rds")
```
