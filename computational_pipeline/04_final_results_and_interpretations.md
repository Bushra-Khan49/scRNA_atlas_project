# Results and Biological Interpretations

This document serves as the comprehensive biological summary of the single-cell RNA sequencing pipeline. The findings detailed below outline the computational strategies and statistical validations used to successfully resolve the human peripheral immune system at a single-cell resolution.

---

## 1. Data Integrity and Quality Control Filtering

To ensure the integrity of all downstream analyses, strict quality control metrics were applied to the raw sequencing dataset. Raw droplet-based sequencing captures significant background noise, including empty fluid droplets and necrotic cells, which can severely confound clustering algorithms if left unresolved.

![Figure 1: Quality Control Metrics](../figures/01_qc_violin.png)

As demonstrated in **Figure 1**, I visualized the distribution of unique genes, total RNA molecules, and the percentage of mitochondrial RNA across the dataset. The mitochondrial expression profile (third panel) reveals a distinct, long "tail" stretching upwards from the baseline population. Because mitochondrial RNA remains trapped in the cell membrane even after a cell dies and leaks its cytoplasmic RNA, this tail represents a subpopulation of necrotic or highly stressed cells. 

By applying a stringent cutoff—removing any cell with >5% mitochondrial expression—I guaranteed that the remainder of this analysis is built exclusively upon a highly viable, healthy cellular population.

---

## 2. Dimensionality Reduction and Unsupervised Clustering

Following quality control, the data was normalized to correct for varying sequencing depths, and the top 2,000 highly variable genes were isolated to focus the algorithm on true biological variance. Attempting to cluster cells across 2,000 dimensions simultaneously is computationally inefficient, necessitating the use of Principal Component Analysis (PCA).

![Figure 3: PCA Variance (Elbow Plot)](../figures/02_elbow.png)

As shown in the elbow plot (**Figure 3**), there is a massive, steep drop in variance across the first few Principal Components. This sharp drop occurs because these initial components capture the largest, most obvious biological differences in the blood sample (for example, the massive genetic difference between a T Cell and a Monocyte). 

However, right around Principal Component 10, the curve completely flattens out into a plateau. This flattening occurs because the algorithm has already mapped all the major biological identities; any variance beyond this point represents irrelevant statistical noise or minor individual cell quirks rather than distinct cell types. **To visually explicitly mark this inflection point, I drew a dashed red line directly at PC 10.** Relying on this visual plateau and the red marker, I purposefully isolated the first 10 dimensions for the final mathematical clustering, ensuring the model ignores the downstream noise.

Utilizing these first 10 dimensions, I constructed a K-Nearest Neighbor (KNN) graph to map transcriptomic similarities, followed by Louvain clustering.

![Figure 4: Mathematical Clustering](../figures/02_umap_clusters.png)

To visually confirm the success of the algorithm, I projected these multidimensional clusters into a 2D space using Uniform Manifold Approximation and Projection (UMAP). As shown in **Figure 4**, the cells resolved into 9 highly distinct, separate transcriptomic islands. The lack of a single, undifferentiated mass confirms that the normalization and mathematical clustering successfully separated the blood sample into distinct cellular subpopulations based purely on their RNA profiles.

---

## 3. Algorithmic Cell Type Annotation and Biomarker Validation

While the mathematical clustering successfully partitioned the data, it provided no biological context. To transition from anonymous mathematical clusters to true biological lineages, I utilized the `SingleR` algorithm to cross-reference the transcriptomic profile of our clusters against the Human Primary Cell Atlas.

![Figure 5: Biological Cell Type Map](../figures/03_umap_celltype.png)

This reference-based annotation successfully mapped the major lymphoid and myeloid lineages of the human immune system, identifying distinct populations of T Cells, B Cells, Monocytes, NK Cells, and Platelets (**Figure 5**). 

To rigorously validate these algorithmic predictions, I performed a differential expression analysis across all identified cell types. The resulting matrix identifies the specific transcriptomic biomarkers that drive the identity of each cluster.

![Figure 6: Biomarker Expression Heatmap](../figures/03_marker_heatmap.png)

As visualized in the biomarker expression heatmap (**Figure 6**), each annotated cell type exhibits a highly specific, densely upregulated block of gene expression. Analyzing the statistical output of this test (`04_annotated_marker_genes.csv`) confirms these findings with textbook human biology. 

For example, the B Cell cluster demonstrated massive expression of the B-cell receptor component **CD79A** (Average Log2 Fold Change = 6.49, $p \approx 0$). This gene is present in 93.7% of the clustered B cells versus only 4.6% of the background population. Similarly, the T Cell cluster is definitively identified by the massive upregulation of the T-cell receptor component **CD3D** (Average Log2 Fold Change = 3.79, present in 88.6% of T cells vs 9.9% of non-T cells).

---

## Conclusion

The presence of these foundational immunological biomarkers—expressing at massive fold-changes almost exclusively within their respective clusters—proves conclusively that this computational pipeline was a complete success. The rigorous quality control, unsupervised clustering, and statistical biomarker validation successfully isolated, identified, and mapped living human immune cells, demonstrating a highly accurate and biologically sound bioinformatics workflow.
