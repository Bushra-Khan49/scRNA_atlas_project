library(Seurat)
library(tidyverse)

# Load the 10x matrix
mat_path <- "data/pbmc3k/filtered_gene_bc_matrices/hg19"
counts   <- Read10X(data.dir = mat_path)

pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k",
                            min.cells = 3, min.features = 200)

# Mitochondrial content = key QC metric (high % = dying/stressed cells)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Premium QC Plots
# Custom colors for the single dot plots to make them pop
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
        pt.size = 0.1, cols = c("#3A86FF", "#8338EC", "#FF006E")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Quality Control Distributions", y = "Count / Percentage")
ggsave("figures/01_qc_violin.png", width = 10, height = 5, bg = "white")

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5, cols = "#3A86FF") +
  labs(title = "Mitochondrial RNA vs Total RNA", 
       x = "Total RNA Molecules Detected (nCount)", 
       y = "Mitochondrial Percentage (%)") +
  theme(legend.position = "none")
ggsave("figures/01_qc_scatter.png", width = 6, height = 5, bg = "white")

# Filtering thresholds — state your reasoning in the report:
# - nFeature_RNA between 200–2500: removes empty droplets and doublets
# - percent.mt < 5%: standard PBMC cutoff, removes dying cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

saveRDS(pbmc, "results/01_pbmc_filtered.rds")
