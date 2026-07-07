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
ggsave("figures/03_marker_heatmap.png", width = 10, height = 8)

# Reference-based automated annotation (industry-standard method, not just manual marker-eyeballing)
ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(pbmc)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

pbmc$cell_type <- pred$labels
DimPlot(pbmc, group.by = "cell_type", label = TRUE)
ggsave("figures/03_umap_celltype.png", width = 7, height = 5)

saveRDS(pbmc, "results/03_pbmc_annotated.rds")
