# CSV File Data Guide

This folder contains the raw tabular data outputs generated during the computational pipeline. This guide explains exactly what these files contain, how to read them, and how they fit into the overall project.

---

## 1. `04_annotated_marker_genes.csv`

### What is in this file?
This file is the master output of the Differential Expression analysis. After we annotated our cell clusters (e.g., T Cells, B Cells, Monocytes), we ran a statistical test to find the "biomarkers"—the specific genes that uniquely identify each cell type. This CSV contains that full list of genes.

### How do you read this table?
When you open this file, you will see several key columns:
- **`gene`**: The official symbol of the gene being tested (e.g., *CD79A*, *LYZ*).
- **`cluster`**: The specific biological cell type that this gene is a marker for (e.g., *B_cell*).
- **`avg_log2FC`**: The "Average Log2 Fold Change". This tells you *how much higher* this gene is expressed in this specific cluster compared to all other cells. A higher number means it is a very strong, highly specific marker.
- **`p_val` & `p_val_adj`**: The statistical significance of the finding. A very low adjusted p-value (e.g., `1e-50`) means we are statistically certain that this gene is a true marker for this cell type, and it's not just random noise.
- **`pct.1` & `pct.2`**: This shows what percentage of cells *in* the cluster (`pct.1`) express the gene, versus the percentage of cells *outside* the cluster (`pct.2`) that express it. A great marker will have a high `pct.1` and a low `pct.2`.

### What does this mean for our work?
This file is the absolute biological proof of our pipeline. While the UMAP images are great for visual presentation, this CSV provides the hard, statistical evidence that our clustering worked. For example, if you look at the `B_cell` rows, you will see classic B cell biology (like CD79A and MS4A1) sitting at the top of the list with massive Fold Changes. This confirms that the algorithm didn't just group random cells together; it successfully isolated true, functional B cells from the blood sample.

### Where does it fit in the pipeline?
This file was generated at the very end of `03_cell_type_annotation.md`. We took the top 5 rows from each cluster in this exact CSV file and used them to draw the final Biomarker Expression Heatmap (Figure 6) to visually showcase the results.
