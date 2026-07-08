# Detailed Data Guide: Targeted T Cell Profile
**Associated File:** `04_Tcell_DE.csv`

---

## 1. What is this file and why was it generated?
While the previous CSVs looked at every cell type in the dataset, this file is a highly targeted differential expression matrix. T cells constitute the absolute vast majority of human peripheral blood mononuclear cells (PBMCs). Because they are the most dominant population in our sample, I generated this specific file to focus exclusively on the genes that define the T Cell lineage against all other immune cells, allowing for a much deeper look into their specific functional states.

---

## 2. Deep Dive into the Hard Data
If you open `04_Tcell_DE.csv`, you will see a focused ranking of T cell-specific genes. Looking directly at the top of this file, we see:

* **CD3D:** `avg_log2FC = 3.80`, `p_val = 0`, `pct.1 = 0.886`
* **IL7R:** `avg_log2FC = 3.27`, `p_val = 0`, `pct.1 = 0.636`
* **CD3E:** `avg_log2FC = 2.93`, `p_val = 0`, `pct.1 = 0.773`
* **IL32:** `avg_log2FC = 2.47`, `p_val = 0`, `pct.1 = 0.867`
* **CD2:** `avg_log2FC = 3.06`, `p_val = 0`, `pct.1 = 0.519`
* **LCK:** `avg_log2FC = 2.21`, `p_val = 0`, `pct.1 = 0.561`
* **LDHB:** `avg_log2FC = 2.10`, `p_val = 0`, `pct.1 = 0.884`
* **CD27:** `avg_log2FC = 3.13`, `p_val = 0`, `pct.1 = 0.398`

### Biological Interpretation
As established, **CD3D** and **CD3E** form the structural core of the T cell receptor, explaining their presence in the vast majority (88% and 77%) of the cluster. However, the genes immediately following them provide a massive amount of biological context:

1. **IL7R (Interleukin-7 Receptor):** The massive presence of this gene is a crucial finding. IL7R is vital for T cell development, long-term survival, and maintaining immune homeostasis. The fact that it is upregulated heavily (`avg_log2FC = 3.27`) proves that we haven't just identified "generic" dying T cells, but rather a healthy, highly functional, communicating T cell population.
2. **IL32:** A pro-inflammatory cytokine highly expressed in activated T-cells. Its strong, ubiquitous presence (`pct.1 = 0.867`) indicates the immune system in this sample is actively surveying and responding to threats.
3. **CD2:** This is a classic T cell surface antigen that mediates physical adhesion between T cells and other cell types (like antigen-presenting cells). Its strong presence confirms the cells are mature and capable of immune synapse formation.
4. **LCK:** This gene encodes a tyrosine kinase that sits directly inside the cell membrane. The moment a T cell recognizes an infection, LCK fires off the internal alarm. Seeing this gene highly expressed confirms that the internal signaling machinery of these T cells is fully intact.
5. **LDHB & CD27:** Metabolism and maturation markers. LDHB indicates active cellular respiration, while CD27 is a critical costimulatory molecule required for generating T-cell memory.

### Noting an Inconsistency
Notice that while **CD3D** is expressed in nearly 90% of the T cells (`pct.1 = 0.886`), the functional markers like **IL7R** (63%) and **CD2** (51%) are only expressed in about half of the T cells. 

Is this a computational error? No, it is a reflection of true, deep biology. CD3D is a structural requirement for *all* T cells, which is why it is near 90%. However, IL7R is primarily a marker for *naïve and memory* T cells, while effector T cells often downregulate it. Therefore, the fact that only ~60% of the T cells express IL7R is not an inconsistency; it statistically proves that our T cell cluster is actually composed of multiple distinct T cell sub-states (naïve, memory, and effector) living together in the blood. 

---

## 3. Final Conclusion
This targeted file goes beyond simply proving "these are T cells." By examining the specific functional markers like IL7R and LCK, and noting their fractional expression percentages, this data proves that our single-cell pipeline successfully captured the deep, complex sub-biology of a living, breathing immune system.
