library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(tidyverse)

# Step 1: Loading the Clustered Data
# I loaded the 2D-mapped, clustered PBMC object from Phase 2. At this point, the clusters only have numbers (0, 1, 2, etc.), and I need to figure out their actual biological identities.
pbmc <- readRDS("results/02_pbmc_clustered.rds")


# Step 2: Marker Gene Discovery (Differential Expression)
# I ran a Wilcoxon Rank Sum test across every single cluster to find genes that are highly expressed in one cluster but turned off in all the others. These are the "biomarkers".
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# I saved this raw mathematical output to a CSV so it can be loaded instantly into the Shiny dashboard.
write.csv(markers, "results/03_marker_genes.csv", row.names = FALSE)


# Step 3: Visualizing the Biomarkers (Heatmap)
# To prove my differential expression worked, I grabbed the top 5 most defining genes for each cluster.
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

# I plotted these in a massive heatmap. I specifically forced the height to 16 inches during the save so that the individual gene names are fully decompressed and legible on the Y-axis.
DoHeatmap(pbmc, features = top5$gene) 
ggsave("figures/03_marker_heatmap.png", width = 12, height = 16, bg = "white")


# Step 4: Automated Biological Annotation
# Instead of manually guessing cell types by eyeballing marker genes (which is prone to human error), I used the industry-standard `SingleR` package to algorithmically match my cells against the Human Primary Cell Atlas database.
ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(pbmc)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)


# Step 5: Final Cell Type Map
# I mapped the algorithm's predictions back onto my original Seurat object.
pbmc$cell_type <- pred$labels

# I generated the final annotated UMAP, again forcing presentation-quality arrows and titles.
DimPlot(pbmc, group.by = "cell_type", label = TRUE, label.size = 4, pt.size = 0.5) +
  theme_classic() +
  labs(title = "Biologically Annotated Cell Types", 
       x = "UMAP Dimension 1", 
       y = "UMAP Dimension 2") +
  theme(
    axis.line = element_line(arrow = arrow(length = unit(0.3, "cm"), type = "closed")),
    plot.title = element_text(face = "bold", size = 14)
  )
ggsave("figures/03_umap_celltype.png", width = 8, height = 6, bg = "white")


# Step 6: Final Data Save
# This fully processed, mathematically clustered, and biologically annotated object is what powers the final interactive dashboard!
saveRDS(pbmc, "results/03_pbmc_annotated.rds")
