library(Seurat)

pbmc <- readRDS("results/03_pbmc_annotated.rds")
Idents(pbmc) <- "cell_type"

# Example: T cells vs everything else
# Since celldex HumanPrimaryCellAtlasData labels T cells as "T_cell" or similar, 
# we check the unique labels first to prevent errors.
cell_types <- unique(pbmc$cell_type)
cat("Available cell types:", paste(cell_types, collapse=", "), "\n")

# Find a T cell related cluster name
tcell_name <- cell_types[grepl("T_cell", cell_types, ignore.case=TRUE)][1]

if(!is.na(tcell_name)) {
    de_results <- FindMarkers(pbmc, ident.1 = tcell_name, ident.2 = NULL,
                               test.use = "wilcox", logfc.threshold = 0.25)
    write.csv(de_results, "results/04_Tcell_DE.csv")
    cat("DE analysis saved to results/04_Tcell_DE.csv\n")
} else {
    cat("Could not find a T cell cluster for DE analysis.\n")
}
