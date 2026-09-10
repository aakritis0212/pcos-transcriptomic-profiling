# PCOS Transcriptomic Profiling & Enrichment Dashboard

An independent computational biology project analyzing differential gene expression and functional pathway enrichment in Polycystic Ovary Syndrome (PCOS). This repository features an interactive web application built with Streamlit and Plotly to explore transcriptomic shifts, Gene Ontology (GO) pathway enrichment, and differential regulation.

## 🧬 Project Overview

Polycystic Ovary Syndrome (PCOS) is a complex metabolic and endocrine disorder. This project processes raw DESeq2 differential expression results to uncover systemic biological trends, highlighting key downregulations in pathways such as cholesterol biosynthesis and mitochondrial translation.

### Key Features
* **Interactive Volcano Plot:** Visualizes genome-wide differential expression (Log2 Fold Change vs. -Log10 Adjusted P-value) with GPU-accelerated rendering (`webgl`) and custom significance mapping.
* **Pre-Ranked GSEA Results:** Evaluates continuous transcriptomic shifts against Gene Ontology biological processes using `gseapy`.
* **Exploreable Data Tables:** Inspect full enrichment outputs, Normalized Enrichment Scores (NES), False Discovery Rates (FDR q-values), and lead genes directly from the browser.

---

## 🛠️ Tech Stack

* **Language:** Python
* **Core Libraries:** `pandas`, `numpy`, `gseapy`, `plotly`
* **Web Framework:** Streamlit (Hosted on Streamlit Cloud)
* **Version Control:** Git & GitHub

---

## 📂 Repository Structure

```text
├── appy.py                           # Main Streamlit dashboard application
├── requirements.txt                  # Python package dependencies
├── pcos_deseq2_results.csv           # DESeq2 differential expression output data
├── Human.GRCh38.p13.annot.tsv        # Human genome annotation reference file
└── gsea_results_precalculated.csv    # Pre-computed Pre-Ranked GSEA results for fast cloud loading
