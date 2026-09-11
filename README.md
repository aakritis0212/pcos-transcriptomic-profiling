🧬 PCOS Transcriptomic Profiling & Regulatory Network Dashboard
An end-to-end bioinformatics pipeline and interactive cloud dashboard for analyzing the transcriptomic drivers of Polycystic Ovary Syndrome (PCOS).

This project transforms raw RNA-seq count matrices from human granulosa cells into actionable biological insights, bridging the gap between raw statistical data and systemic biological understanding.

🔗 Launch the Live Interactive Dashboard Here

💡 Why You Should Explore This Project
For Biological Researchers:
PCOS is a highly complex endocrine and metabolic disorder. Rather than just looking at isolated gene changes, this pipeline evaluates the entire system. It identifies exactly which biological pathways are malfunctioning (via Gene Ontology), flags the specific transcription factors orchestrating these changes (via ChIP-X screening), and maps out physical protein interactions to identify potential therapeutic hubs.

For Data Scientists & Bioinformaticians:
This repository demonstrates a robust, production-ready data science workflow that handles real-world biological data complexities:

Strict Data Validation: Automatically detects and resolves file encoding issues (UTF-16) and filters out multi-assay contamination (removing miRNA/lncRNA noise) to guarantee a mathematically pristine mRNA matrix.

Pre-computation Architecture: Heavy statistical processing (PyDESeq2, 1000-permutation GSEA) is executed locally. The cloud dashboard relies on ultra-lightweight precalculated CSVs, resulting in a blazing-fast user experience.

Dynamic API Integration: The app dynamically queries the STRING database API to build real-time network graphs without storing massive interaction databases locally.

📊 What the Results Signify
The dashboard is divided into four main analytical tiers. Here is what the data actually tells us:

1. Differential Expression (The "What")
Visuals: Volcano Plot, MA Plot, Z-Score Heatmap.

Significance: By utilizing the DESeq2 algorithm on strictly validated mRNA samples (5 PCOS vs 5 Healthy Controls), we calculate the exact Log2 Fold Change of every gene. This reveals which specific genes are aggressively up-regulated (turned on) or down-regulated (turned off) in the granulosa cells of PCOS patients.

2. Functional Pathways via GSEA (The "How")
Visuals: Enrichment Bar Plots, Lead-Gene Bubble Plots.

Significance: Genes do not act alone. Gene Set Enrichment Analysis (GSEA) evaluates our ranked list against Gene Ontology (GO) Biological Processes. Instead of saying "Gene X is up," this tells us "The entire pathway for lipid metabolism is disrupted" or "Inflammatory cascades are hyperactive," providing the actual biological narrative of the disease.

3. Master Regulators / TFs (The "Who")
Visuals: ChIP-X Transcription Factor Bar Plot.

Significance: This is the most crucial regulatory step. By screening our data against the ENCODE ChIP-X database, we identify the upstream Transcription Factors (TFs). These are the "puppeteers" driving the genetic changes. If we want to understand what causes the widespread dysregulation in PCOS, these TFs are the prime suspects.

4. Protein-Protein Interaction (PPI) Network (The "Structure")
Visuals: Interactive NetworkX / Plotly Node Graph.

Significance: We extract the "core driver genes" from our top disrupted pathways and query the STRING database. The resulting network shows how these proteins physically interact in 3D space. Highly connected nodes (large hubs) represent critical structural pinch-points in the cell—often making them the best targets for pharmacological drugs.

🛠️ Data Origin & Methodology
Dataset: NCBI GEO GSE168404.

Filtering Protocol: The original repository contained multi-assay data. This pipeline enforces a strict 10-sample filter to isolate purely mRNA assays:

Controls (n=5): C1, C2, C3, C4, C5

PCOS (n=5): P1, P2, P3, P4, P5

Differential Expression: PyDESeq2

Enrichment: GseaPy (Pre-ranked GSEA, 1000 permutations)

Visualization: Plotly Express, NetworkX, Streamlit


├── GSE168404_mRNA_C-vs-P.all.txt   # Raw, UTF-16 validated mRNA count matrix
├── Human.GRCh38.p13.annot.tsv      # Gene ID to HUGO symbol mapping reference
│
# --- PHASE 1: LOCAL STATISTICAL PIPELINE (Do not run on cloud) ---
├── differentialexpression.py       # Cleans data, runs PyDESeq2, outputs DESeq2 stats
├── enrichmentanalysis.py           # Runs GSEA for GO Pathways and TF ChIP-X screening
│
# --- PHASE 2: STREAMLIT DASHBOARD DEPLOYMENT ---
├── appy.py                         # The main Streamlit dashboard application
├── pcos_deseq2_results.csv         # PRE-CALCULATED: DESeq2 output
├── gsea_results_precalculated.csv  # PRE-CALCULATED: Pathway enrichment output
├── tf_screening_results.csv        # PRE-CALCULATED: Transcription factor output
├── requirements.txt                # Cloud deployment dependencies


