# Single-Cell RNA-Seq Transcriptomics Explorer (PBMC) 🧬

![Shiny App Demo](https://raw.githubusercontent.com/plotly/plotly.js/master/logo/plotlyjs-logo.svg) 
*(Note: Replace with a screenshot of the actual dashboard!)*

## 🚀 Why I Built This
I built this as a solo, deep-dive individual project to pressure-test my R programming skills and prove my ability to handle end-to-end Computational Biology workflows. I wanted to go beyond just running standard tutorials—I wanted to take raw, noisy 10x Genomics sequencing data and manually sculpt it into a stunning, production-ready SaaS dashboard that anyone can explore. 

This repository documents my entire analytical journey: from the messy early stages of data wrangling, to the complex mathematics of dimensionality reduction, all the way to deploying a responsive Shiny UI.

---

## 🔬 The Step-by-Step Analytical Walkthrough

Here is exactly how I processed the data, and the biological reasoning behind every single decision I made.

### Phase 1: Quality Control (QC) & Data Scrubbing
You can't do good science with garbage data. When sequencing cells, some cells die, rupture, or get trapped together (doublets). I needed to clean this up mathematically.

- **The Mitochondrial Cutoff (`< 5%`):** When a cell ruptures, its cytoplasmic RNA leaks out, but the RNA trapped inside its mitochondria remains. Therefore, a massive spike in mitochondrial RNA means the cell is dying. I used a regular expression (`^MT-`) to calculate the mitochondrial percentage of every cell and aggressively filtered out anything above 5%.
- **The Gene Count Threshold (200 - 2500):** If a cell had less than 200 genes, it was likely an empty droplet of fluid. If it had more than 2,500, it was probably two cells stuck together. I trashed both.

![QC Metrics](figures/01_qc_violin.png)

### Phase 2: Dimensionality Reduction (PCA & UMAP)
Human blood cells express roughly 20,000 genes. Trying to visualize a 20,000-dimensional matrix is impossible. 

- **PCA:** First, I ran Principal Component Analysis to compress those 20,000 dimensions down to the top 10 vectors of highest mathematical variance.
- **UMAP:** Then, I used UMAP (Uniform Manifold Approximation and Projection) to crush those 10 dimensions down to a 2D map. This allows us to visually cluster cells that share similar biological "fingerprints." 

![UMAP Clusters](figures/02_umap_clusters.png)

### Phase 3: Differential Expression & Biomarker Discovery
Once the cells were clustered in 2D space, I needed to figure out what they actually were. I ran a **Wilcoxon Rank Sum test** (`FindAllMarkers`) to find the statistically significant genes defining each cluster.

> **💡 Fun Biology Fact:** During this phase, I noticed genes like *ATG12* (an autophagy gene) were completely missing from my marker tables, despite being highly expressed in the raw data. Why? Because *ATG12* is a basic survival "housekeeping" gene. Every cell expresses it. Because it doesn't uniquely define a *single* cell type, my algorithm correctly threw it out!

Finally, I used reference-based automated annotation (`SingleR` + `celldex`) to map these clusters to true human immune identities (T-cells, B-cells, etc.).

![Annotated Cell Types](figures/03_umap_celltype.png)
![Marker Heatmap](figures/03_marker_heatmap.png)

### Phase 4: The Interactive Dashboard
I didn't just want a static report; I wanted a fully interactive explorer. I built a premium R Shiny app utilizing `bslib` for a stunning Light/Dark mode toggle and `plotly` for WebGL-powered interactive charts that render thousands of cells without lagging the browser. 

---

## 🛠️ Challenges & Troubleshooting
This project definitely threw some curveballs at me. Here is how I solved them:

1. **The C++ Compilation Nightmare (Phase 1):** While installing Seurat and `hdf5r`, I ran into brutal Mac-specific C++ compiler errors because R was trying to build the packages from source without the correct Fortran/C binaries. I had to step out of R, configure my system's Homebrew compiler pathways, and force binary installations to finally get the environment stabilized.
2. **Bioconductor Dependency Hell (Phase 3):** Getting `SingleR` and `celldex` to talk to each other was tricky. Bioconductor has incredibly strict versioning rules, and mismatched dependencies were causing silent failures when querying the Human Primary Cell Atlas. I ended up completely wiping my package cache and orchestrating a fresh installation via `BiocManager` to ensure version harmony.
3. **The 'FetchData' UI Crash (Phase 4):** When building the Shiny dashboard, I noticed that if a user rapidly deleted the gene name from the search box, the app would instantly crash. The server was trying to search the Seurat object for a blank string (`""`). I patched this by injecting a strict validation logic (`req(input$gene %in% rownames(pbmc))`) to ensure the server only attempts to render physically existing genes.
4. **GitHub's 100MB File Limit:** Single-cell datasets are massive. My finalized `.rds` files were nearly 300MB! When I tried to push this to GitHub, the terminal locked up. I had to manually intercept the Git index, untrack the massive matrices using `git rm --cached`, and rewrite my `.gitignore` to protect the repository from crashing while still preserving my actual code.

---

## 💻 Try it Yourself (Reproducibility)

I locked all 231 dependencies using `renv` so this project is 100% reproducible. 

To run this on your own machine:
1. Clone the repo: `git clone https://github.com/Bushra-Khan49/scRNA_atlas_project.git`
2. Open it in RStudio.
3. Restore the environment: `renv::restore()`
4. Run the app: `shiny::runApp("shiny_app")`
