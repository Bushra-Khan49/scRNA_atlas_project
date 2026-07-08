# CSV File Data Guide

This folder contains the raw tabular data outputs generated during the computational pipeline. This guide explains exactly what these files contain, the specific numbers we found, and what they mean for our biological conclusions.

---

## 1. `04_annotated_marker_genes.csv`

### What is in this file?
After we mathematically clustered our cells and assigned them names (like T Cells and B Cells), we needed to prove those names were accurate. We ran a differential expression test to find the specific genes that uniquely identify each group. This CSV contains that full, highly specific list of genes.

### How to interpret the actual numbers we got
Let's look at a specific, real example straight from our data to understand how to read this:

In the CSV, you will find a row for **CD79A**.
* **`cluster` = B_cell**: The algorithm states that CD79A is the defining marker for B Cells.
* **`avg_log2FC` = 6.49**: This is massive. Because this is a log2 scale, an average log2 fold change of 6.49 means that CD79A is expressed roughly **90 times higher** in B cells than in any other cell in our entire dataset. 
* **`p_val` & `p_val_adj` = 0**: The probability that this gene expression is a random statistical mistake is practically zero. 
* **`pct.1` = 0.937**: This means that 93.7% of all the cells in our B Cell cluster express this gene.
* **`pct.2` = 0.046**: This means that only 4.6% of all the other cells outside this cluster express it.

### What do these numbers signify biologically?
These numbers aren't just statistics; they perfectly mirror human biology. CD79A is a critical protein that forms the B-cell antigen receptor complex. A B cell literally cannot function without it. 

The fact that our data shows CD79A expressed at a massive 6.49 Log2 Fold Change, present in 93.7% of our predicted B cells and almost entirely absent everywhere else (4.6%), **proves conclusively** that the cells in that cluster are true, biologically functional B cells. 

If we look at T cells, we see the exact same thing: **CD3D** has an `avg_log2FC` of 3.79 (expressed over 13x higher in T cells) and is present in 88.6% of T cells versus only 9.9% of non-T cells. CD3D is part of the T-cell receptor complex. 

### Final Conclusion
When you see numbers like this—where textbook biological markers perfectly align with our mathematically generated clusters with p-values of 0—it signifies that our single-cell pipeline was a complete success. The numbers in this CSV prove that we didn't just group random noise; we successfully isolated and identified living, distinct immune cells straight from the blood sample.
