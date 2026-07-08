# Single-Cell RNA-Seq Transcriptomics Explorer (PBMC)

This repository contains my full end-to-end Computational Biology pipeline and an interactive visualization dashboard built in R. I process raw Peripheral Blood Mononuclear Cell (PBMC) scRNA-seq data through quality control, normalization, dimensionality reduction, cell-type annotation, and deploy it as a premium interactive SaaS-style Shiny application.

## 🧬 Biological Rationale & Methods

### 1. Quality Control (QC)
In single-cell RNA sequencing, it is vital to mathematically remove dead, dying, or ruptured cells before analysis. 
- **Mitochondrial Threshold (< 5%):** Apoptotic or lysed cells leak their cytoplasmic RNA, leaving behind a disproportionate ratio of mitochondrial RNA (which is protected inside the mitochondrial membrane). I strictly filtered out cells with `> 5%` mitochondrial RNA.
- **Gene Count Threshold (200 < nFeature_RNA < 2500):** I excluded empty droplets (cells expressing very few genes) and doublets (two cells captured in one droplet, resulting in abnormally high gene counts).

### 2. Dimensionality Reduction (PCA & UMAP)
Human cells express ~20,000 genes, representing 20,000 mathematical dimensions. 
- **PCA (Principal Component Analysis):** I reduced this to the top 10 most statistically variable principal components to capture the highest biological variance.
- **UMAP (Uniform Manifold Approximation and Projection):** I compressed the 10 dimensions down to 2D space to visually cluster biologically identical immune cells.

### 3. Differential Expression & Biomarker Annotation
I utilized the Wilcoxon Rank Sum test via Seurat's `FindAllMarkers` algorithm to identify statistically significant biomarker genes that define each cluster. 
- *Note on Housekeeping Genes:* Basic survival genes (e.g., *ATG12*, an autophagy-related gene) are expressed across all cells. They are excluded from my Biomarker tables because they do not uniquely separate distinct immune identities, proving the pipeline's statistical rigor.
- Clusters were then biologically annotated into true cell identities (e.g., T-cells, B-cells, Monocytes) using reference mapping.

## 🚀 The Interactive Dashboard

The results are compiled into a premium R Shiny Dashboard built with `bslib` and `plotly`.

### Features:
- **Liquid Premium UI:** Fully responsive dashboard featuring dynamic Light/Dark mode tailored to high-end SaaS standards.
- **Server-Side Selectize:** Employs `server = TRUE` to instantly load all 13,000+ sequenced genes without lagging the browser.
- **Dynamic Plotly WebGL:** Renders thousands of individual single cells effortlessly with custom color palettes and interactive hover states.

## 🛠️ Reproducible Environment (`renv`)

This project is fully production-ready and reproducible. All dependencies (Seurat, Shiny, Bslib, Plotly, DT) are tracked via `renv`.

### How to Run Locally

1. Clone this repository.
2. Open the project in RStudio.
3. Restore the exact package environment:
   ```R
   renv::restore()
   ```
4. Launch the dashboard:
   ```R
   shiny::runApp()
   ```

## 📂 Project Architecture

```
scRNA_atlas_project/
├── renv/                       # Isolated reproducible R environment
├── renv.lock                   # Package versions and hashes
├── data/                       # Raw cellranger output (matrix, features, barcodes)
├── computational_pipeline/     # My detailed step-by-step pipeline tutorials (.md)
├── results/                    # Compiled .rds objects and CSV biomarker data
├── figures/                    # Custom generated ggplot2 visualizations
├── www/                        # Custom CSS (custom.css)
├── app.R                       # Full Dashboard architecture script
└── README.md                   # Documentation
```
