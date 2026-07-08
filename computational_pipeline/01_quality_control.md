# Project Context: Mapping the Human Immune System

Before diving into the code, it is critical to understand the overarching goal of this project. The human peripheral immune system is incredibly complex, composed of many different types of cells (T cells, B cells, monocytes, etc.) all circulating together in the blood. Our objective is to take a raw blood sample, sequence the RNA of every single cell individually using 10x Genomics technology, and computationally reconstruct a high-resolution map of the immune system. 

By analyzing the unique transcriptomic signature of each cell, we can computationally group them, identify them, and validate them. This file, **Phase 1**, is the foundation of that entire process. Before we can map the cells, we must mathematically separate the living, healthy cells from the dead cells and sequencing errors.

---

# Phase 1: Quality Control

In this initial phase, I process the raw single-cell RNA sequencing data. Raw transcriptomic data contains significant background noise, including empty fluid droplets (where the machine captured no cell) and necrotic cells (dead cells). If we don't filter these out immediately, they will ruin the downstream clustering.

## 1. Environment Setup

It is important to manage dependencies carefully when working with single-cell datasets. Here, I initialize the environment by loading `Seurat`, which serves as the core computational framework for this analysis, alongside `tidyverse` for data manipulation, and `patchwork` to perfectly stitch our visualizations together later.

```r
# Load the Seurat package for single-cell analysis
library(Seurat)
# Load tidyverse for data wrangling
library(tidyverse)
# Load ggplot2 and patchwork for custom visualizations
library(ggplot2)
library(patchwork)
```

---

## 2. Loading the Raw Data

The sequencing data is provided as a standard 10x Genomics output (comprising matrix, barcodes, and features files). I used the `Read10X` function to parse these files into a sparse matrix, which is then initialized into a `SeuratObject`. 

During initialization, I applied a lenient baseline filter (`min.cells = 3` and `min.features = 200`). This is a critical best practice to immediately discard unexpressed genes and nearly empty droplets, which significantly optimizes memory usage for the rest of the pipeline.

```r
# Point to the directory containing the raw matrix files
mat_path <- "../data/pbmc3k/filtered_gene_bc_matrices/hg19"

# Load the raw 10x data into a sparse matrix
counts   <- Read10X(data.dir = mat_path)

# Create the Seurat object and apply a lenient initial filter
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)
```

---

## 3. Quantifying Mitochondrial RNA

A primary indicator of cell death in single-cell data is an abnormally high ratio of mitochondrial RNA. When a cell dies and its membrane ruptures, its normal cytoplasmic RNA leaks out, but the tough mitochondria remain trapped inside. 

To identify these compromised cells, I calculated the percentage of total RNA counts originating from mitochondrial genes (identified by the `MT-` prefix).

```r
# Calculate the percentage of counts originating from mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

---

## 4. Visualizing Quality Metrics

Before establishing filtering thresholds, I must visualize the distribution of the data to understand the baseline quality of the sample. Note that the code below is exactly how I built the custom 3-panel figure, explicitly adding arrows to the axes and placing the title at the bottom.

```r
# Generate the first plot: Unique Genes
p1 <- VlnPlot(pbmc, features = "nFeature_RNA", pt.size = 0.1, cols = "#3A86FF") + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"), type = "closed")), legend.position="none", plot.title=element_text(size=12, face="bold", hjust=0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(title="Unique Genes Detected", x="", y="Number of Genes")

# Generate the second plot: Total RNA Count
p2 <- VlnPlot(pbmc, features = "nCount_RNA", pt.size = 0.1, cols = "#FF006E") + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"), type = "closed")), legend.position="none", plot.title=element_text(size=12, face="bold", hjust=0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(title="Total RNA Molecules", x="", y="RNA Count")

# Generate the third plot: Mitochondrial Percentage
p3 <- VlnPlot(pbmc, features = "percent.mt", pt.size = 0.1, cols = "#8338EC") + 
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"), type = "closed")), legend.position="none", plot.title=element_text(size=12, face="bold", hjust=0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(title="Mitochondrial RNA", x="", y="Percentage (%)")

# Stitch them together using patchwork and place the overall title at the bottom
final_vln <- (p1 | p2 | p3) + plot_annotation(caption = "Figure 1: Quality Control Metrics (Pre-Filter)", theme = theme(plot.caption = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(t = 20))))

# Save the custom image
ggsave("../figures/01_qc_violin.png", plot = final_vln, width = 12, height = 5, bg = "white")
```

![Figure 1: Quality Control Metrics](../figures/01_qc_violin.png)

### Image Interpretation: How I inferred the cell quality
If you look closely at **Figure 1**, specifically the third panel labeled "Mitochondrial RNA", you can see how I made my filtering decisions purely based on the physical shape of the graph. 

Because the width of a violin plot represents the raw *number* of cells, you can see a massive, dense bulge at the very bottom between 0% and 5%. This massive width visually proves that the vast majority of the cells in our sample are perfectly healthy with low mitochondrial leakage. However, if you look directly above the 5% mark, you see a very thin "tail" stretching upwards. Because that tail is so thin compared to the bottom bulge, I can definitively infer that these highly stressed/dead cells make up a very small minority of the total population. 

---

## 5. Applying Biological Cutoffs

Relying entirely on the visual evidence from the violin plot above, I applied a strict filtering threshold to cut off that long tail (`percent.mt < 5`). 

Additionally, looking at the first panel ("Unique Genes Detected"), I can see a thin tail stretching above 2,500 genes. These are likely doublets (two cells accidentally captured in one droplet), so I removed them. I also set a lower bound of 200 genes to trim the very bottom tip of empty droplets.

```r
# Filter the Seurat object to keep only high-quality, viable cells based on the visual cutoffs
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

---

## 6. Saving the Processed Dataset

Because processing large genomic matrices is computationally expensive, I saved the filtered Seurat object as an `.rds` file. This checkpoint allows me to directly load the clean dataset into Phase 2 without repeating these quality control steps.

```r
# Save the filtered Seurat object as an R data file
saveRDS(pbmc, "../results/01_pbmc_filtered.rds")
```

---

### Phase 1 Complete
**Conclusion & Next Steps:** We have successfully loaded the raw dataset and applied strict biological cutoffs to filter out necrotic cells and empty droplets. With a clean, highly viable matrix of human immune cells now saved to disk, we are ready to proceed to the next stage of the pipeline. In **Phase 2: Normalization and Clustering**, we will take these healthy cells, normalize their sequencing depths, and mathematically group them into distinct transcriptomic islands.
