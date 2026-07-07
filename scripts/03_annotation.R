library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)

pbmc <- readRDS("results/02_pbmc_clustered.rds")

# Marker gene discovery per cluster (Wilcoxon test, same statistical family as ROC/PR-AUC based tools)
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "results/03_marker_genes.csv", row.names = FALSE)

top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
DoHeatmap(pbmc, features = top5$gene) 
ggsave("figures/03_marker_heatmap.png", width = 12, height = 16, bg = "white")

# Reference-based automated annotation (industry-standard method, not just manual marker-eyeballing)
ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(pbmc)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

pbmc$cell_type <- pred$labels
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

saveRDS(pbmc, "results/03_pbmc_annotated.rds")
