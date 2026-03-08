<<<<<<< HEAD
# RNA-Seq Differential Expression Pipeline + Shiny Dashboard

![R](https://img.shields.io/badge/R-4.5.2-276DC3?style=flat&logo=r)
![DESeq2](https://img.shields.io/badge/DESeq2-Bioconductor-brightgreen)
![Shiny](https://img.shields.io/badge/Shiny-Dashboard-blue)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Dataset](https://img.shields.io/badge/GEO-GSE157234-red)

A complete, reproducible RNA-Seq differential expression analysis pipeline with an interactive Shiny dashboard. Built using real published data from **Shemer et al., Immunity 2020**.

---

## 🚀 Live Demo

**[Launch Dashboard →](https://connect.posit.cloud/rna-seq-shiny-pipeline)**

> ⚠️ First load takes ~2 minutes — DESeq2 is running live on real data!

---

## 🧬 Biological Question

> What happens to microglia when they lose the ability to sense IL-10 after an immune challenge?

| | |
|---|---|
| **Dataset** | GSE157234 — Mouse microglia, 48h after peripheral LPS challenge |
| **Comparison** | IL-10R Mutant (deficient) vs Control (proficient) microglia |
| **Key finding** | Without IL-10 signalling, microglia hyperactivate and overproduce TNF, causing neuronal damage |
| **Paper** | Shemer et al., *Immunity* 53, 1033–1049, 2020 |

---

## 📊 Dashboard Features

| Feature | Description |
|---|---|
| 🌋 **Volcano Plot** | Interactive — hover any gene, adjust padj and LFC thresholds live |
| 🔵 **PCA Plot** | Sample clustering — confirms Mutant vs Control separation at 48h |
| 🟥 **Heatmap** | Top N DEGs with z-scored expression, adjustable gene count |
| 📋 **Results Table** | Searchable, filterable DEG table with CSV download |
| 📤 **Upload Your Data** | Upload your own count matrix + metadata to reuse the full pipeline |
| ⬇️ **Downloads** | PNG, PDF, and CSV exports for all plots and results |

---

## 🖼️ Screenshots

### Volcano Plot
![Volcano Plot](plots/volcano_plot.png)

### PCA Plot
![PCA Plot](plots/pca_plot.png)

### Heatmap
![Heatmap](plots/heatmap_top50_DEGs.png)

---

## 📁 Repository Structure

```
rna-seq-shiny-pipeline/
│
├── README.md
├── .gitignore
│
├── files/                                  # Local analysis (run in RStudio)
│   ├── analysis_final.R                    # Full DESeq2 pipeline
│   ├── app_final.R                         # Shiny app (local version)
│   ├── data/
│   │   ├── count_matrix_48h_clean.csv      # Processed count matrix (48h)
│   │   └── metadata_48h_clean.csv          # Sample metadata
│   ├── results/
│   │   ├── DESeq2_results_Mutant_vs_Control.csv
│   │   ├── top100_upregulated.csv
│   │   └── top100_downregulated.csv
│   └── plots/
│       ├── volcano_plot.png
│       ├── pca_plot.png
│       └── heatmap_top50_DEGs.png
│
└── deploy/                                 # Deployment version (Posit Cloud)
    ├── app.R                               # Shiny app (deployment version)
    ├── manifest.json                       # Package versions for Posit Cloud
    └── data/
        ├── count_matrix_48h_clean.csv
        └── metadata_48h_clean.csv
```

---

## ⚙️ How to Run Locally

### 1. Clone the repository
```bash
git clone https://github.com/mdabrarfaiyaj/rna-seq-shiny-pipeline.git
cd rna-seq-shiny-pipeline
```

### 2. Install required packages
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery", "SummarizedExperiment"))

install.packages(c("shiny", "shinydashboard", "ggplot2", "ggrepel",
                   "pheatmap", "dplyr", "RColorBrewer", "plotly", "DT"))
```

### 3. Run the analysis pipeline
```r
# Downloads GEO data, runs DESeq2, saves all result files
source("files/analysis_final.R")
```

### 4. Launch the Shiny dashboard
```r
shiny::runApp("files/app_final.R", launch.browser = TRUE)
```

---

## 🔬 Methods

| Step | Tool | Details |
|---|---|---|
| Data download | GEOquery | GSE157234 supplementary files |
| Input data | UTAP pipeline | Normalized counts, rounded to integers for DESeq2 |
| Sample subset | Manual | 48h post-LPS only — peak of hyperactivation |
| Excluded samples | Manual | DKO (double-knockout) — third genotype, not part of comparison |
| Differential expression | DESeq2 | design = `~ condition` |
| Low-count filtering | DESeq2 | ≥10 counts in ≥2 samples |
| Significance | DESeq2 | padj < 0.05 AND \|log2FC\| > 1 |
| Transformation | DESeq2 VST | Variance-stabilizing transform for PCA and heatmap |
| Visualisation | ggplot2, pheatmap, plotly | Volcano, PCA, Heatmap |

---

## 📈 Key Results

| Direction | Genes | Biological Meaning |
|---|---|---|
| ⬆️ Upregulated in Mutant | *Tnf*, *Ccl5*, *Il12b*, *Il6*, *Il1b* | Pro-inflammatory hyperactivation |
| ⬇️ Downregulated in Mutant | *P2ry12*, *Sall1*, *Tmem119* | Loss of homeostatic microglia identity |

**Conclusion:** Loss of IL-10 receptor signalling prevents microglia from returning to ground state after LPS challenge — consistent with Figure 3 of the paper (~954 upregulated, ~693 downregulated genes at 48h).

---

## 📤 Use This Pipeline for Your Own Data

The Shiny dashboard accepts custom uploads:
- **Count matrix** — CSV, rows = genes, columns = samples, raw integer counts
- **Metadata** — CSV, rows = samples, must include a `condition` column with exactly 2 groups

---

## 👤 Author

**Md. Abrar Faiyaj** — Bioinformatics Analyst

[![GitHub](https://img.shields.io/badge/GitHub-mdabrarfaiyaj-black?logo=github)](https://github.com/mdabrarfaiyaj)

🔬 **Services:** Bulk RNA-Seq (DESeq2, edgeR) · scRNA-Seq (Seurat) · Custom Shiny Dashboards · Pathway Enrichment (GSEA, GO, KEGG)

---

## 📄 Dataset Reference

**GEO:** [GSE157234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157234)  
**Paper:** Shemer A et al. *Interleukin-10 Prevents Pathological Microglia Hyperactivation following Peripheral Endotoxin Challenge.* Immunity. 2020;53(5):1033–1049.  
**DOI:** [10.1016/j.immuni.2020.09.018](https://doi.org/10.1016/j.immuni.2020.09.018)

---

## 📜 License

MIT License — free to use and adapt with attribution.

