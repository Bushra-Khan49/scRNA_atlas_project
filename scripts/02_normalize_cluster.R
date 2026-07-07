library(Seurat)
library(ggplot2)

pbmc <- readRDS("results/01_pbmc_filtered.rds")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

ElbowPlot(pbmc, ndims = 30)
ggsave("figures/02_elbow.png", width = 5, height = 4)
# Pick number of PCs where the "elbow" flattens — document your choice (commonly 10-15 for PBMC3k)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)   # Louvain algorithm
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  theme_classic() +
  labs(title = "Mathematical Clustering (UMAP)", 
       x = "UMAP Dimension 1", 
       y = "UMAP Dimension 2") +
  theme(
    axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
    plot.title = element_text(face = "bold", size = 14)
  )
ggsave("figures/02_umap_clusters.png", width = 7, height = 6, bg = "white")

saveRDS(pbmc, "results/02_pbmc_clustered.rds")
