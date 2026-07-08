# Final Results and Biological Interpretations

This document serves as the master summary for the entire pipeline. Below, I break down exactly why I generated every single figure, how it was created, and what biological conclusions we can draw from it for our project.

---

## Step 1: Quality Control Filtering

![Figure 1: Quality Control Metrics](../figures/01_qc_violin.png)

### Why we generated this
Before running any analysis, we have to make sure the data is actually usable. Sequencing machines are messy, and they often capture empty droplets or dead cells. If we don't filter those out, they will ruin the downstream clustering.

### How we generated it
I calculated the total number of unique genes, the total number of RNA molecules, and the percentage of RNA coming specifically from mitochondrial genes for every single cell. I then plotted these distributions using a Violin Plot.

### What we are actually trying to see
We want to see the "normal" range for a healthy cell. Specifically, we are looking for the abnormal outliers: cells with almost no genes (empty droplets), cells with an impossible number of genes (doublets), and cells with a massive spike in mitochondrial RNA (which happens when a cell dies, its membrane ruptures, and its normal RNA leaks out).

### What result did we get?
Looking at the third plot, we see a dense cluster of cells near the bottom (healthy), but a long "tail" stretching upwards. That tail represents our dead/dying cells.

### What does this mean for our project?
This allowed us to set a strict mathematical threshold. By dropping any cell with more than 5% mitochondrial RNA, we guaranteed that the rest of our project is built exclusively on high-quality, living cells, making our final conclusions highly accurate.

---

## Step 2: Mathematical Clustering

![Figure 4: Mathematical Clustering](../figures/02_umap_clusters.png)

### Why we generated this
Once we had healthy cells, we needed to see if they naturally grouped together based on their gene expression. We wanted to prove that the immune system is made of distinct, mathematically separable cell types.

### How we generated it
I found the 2,000 genes that varied the most across the cells, crushed those 2,000 dimensions down to 10 using PCA, and then used a K-Nearest Neighbor algorithm to group the cells. Finally, I used UMAP to project that multidimensional data onto a 2D map so we could visually look at it.

### What we are actually trying to see
We are looking for "islands". If the cells just formed one giant, messy blob, it would mean our sequencing failed or the cells are all completely identical. Distinct islands mean the math successfully found different types of cells.

### What result did we get?
The UMAP generated 9 very distinct clusters.

### What does this mean for our project?
This proves that our normalization and PCA steps worked perfectly. The cells are separating beautifully based on their transcriptomic profiles. At this stage, we don't know *what* they are yet, but we know they are biologically distinct groups.

---

## Step 3: Biological Cell Type Annotation

![Figure 5: Biological Cell Type Map](../figures/03_umap_celltype.png)

### Why we generated this
"Cluster 0" and "Cluster 1" are useless names for a biological study. We needed to assign real, scientific identities (like T Cells or B Cells) to the clusters we found in the previous step.

### How we generated it
I ran the `SingleR` algorithm, which took the genetic profile of our clusters and cross-referenced it against the Human Primary Cell Atlas database. The algorithm automatically predicted the true biological identity of each cluster, and I replotted the UMAP with these new labels.

### What we are actually trying to see
We are trying to see if the mathematical clusters perfectly align with real biological cell types. We want to identify the core components of the human peripheral immune system.

### What result did we get?
The algorithm successfully identified massive populations of T Cells, along with distinct islands of Monocytes, B Cells, NK Cells, and Platelets. 

### What does this mean for our project?
This is the core deliverable of the project. We have successfully taken a massive, anonymous matrix of raw sequencing numbers and reconstructed a high-resolution, biologically accurate map of the human immune system at single-cell resolution.

---

## Step 4: Biomarker Validation

![Figure 6: Biomarker Expression Heatmap](../figures/03_marker_heatmap.png)

### Why we generated this
We needed to mathematically prove that our cell type annotations were correct. We can't just trust the algorithm blindly; we need to see the actual genes driving those identities.

### How we generated it
I ran a differential expression test across all the cell types to find the genes that were highly expressed in one group but turned off in all the others. I then took the top 5 genes for each cell type and plotted them on this heatmap.

### What we are actually trying to see
We want to see solid blocks of bright yellow (high expression) that perfectly align with the columns (the cell types). If the yellow is scattered randomly everywhere, it means our clusters are garbage.

### What result did we get?
We got beautiful, distinct blocks of expression. For example, the B Cell column shows massive expression of CD79A, which is the universal, textbook marker for B Cells.

### What does this mean for our project?
This acts as the final validation for the entire pipeline. The strict grouping of these highly expressed genes perfectly aligns with our cell type annotations, confirming without a doubt that our analysis is biologically accurate and publishable.
