library(Seurat)
library(ggplot2)

# Step 1: Loading the Clean Data
# I loaded the filtered PBMC dataset from the QC phase so I could begin the mathematical transformations.
pbmc <- readRDS("results/01_pbmc_filtered.rds")

# Step 2: Normalization & Feature Selection
# Raw RNA counts are highly skewed by sequencing depth (some cells get sequenced deeper than others).
# I normalized the data using LogNormalize (scaling by 10,000 to match standard transcripts-per-million rates).
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Out of ~20,000 genes, most are just background noise. 
# I ran a variance stabilizing transformation (VST) to isolate the top 2000 most highly variable genes 
# that actually drive biological differences between cell types.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# I centered and scaled the data so highly expressed genes don't dominate the downstream PCA math.
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Step 3: Dimensionality Reduction (PCA)
# A 2000-dimension matrix is computationally impossible to cluster effectively.
# I ran Principal Component Analysis (PCA) to compress the data down to linear components of maximum variance.
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# I generated an Elbow Plot to visualize the drop-off in variance. 
# I visually inspected this and determined that the vast majority of biological signal was captured in the first 10 PCs.
ElbowPlot(pbmc, ndims = 30)
ggsave("figures/02_elbow.png", width = 5, height = 4, bg = "white")

# Step 4: Manifold Approximation & Clustering
# Using those top 10 PCs, I built a K-Nearest Neighbor (KNN) graph to map cell-to-cell proximities.
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# I applied the Louvain algorithm (resolution = 0.5) to group these neighbors into distinct biological clusters.
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Finally, I crushed the 10-dimensional space down to 2 dimensions using UMAP so I could physically see the clusters.
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Step 5: Visualizing the Clusters
# I generated the UMAP plot. I specifically injected theme_classic() to override Seurat's default 
# boxed aesthetic, applying professional X/Y axis arrows and bold titles for presentation-quality output.
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

# Step 6: Save Checkpoint
# I saved the fully clustered object to pass to the final annotation phase.
saveRDS(pbmc, "results/02_pbmc_clustered.rds")
