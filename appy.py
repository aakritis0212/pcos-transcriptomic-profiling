import streamlit as st
import pandas as pd
import gseapy as gp
import plotly.express as px

st.set_page_config(page_title="PCOS Transcriptomics", layout="wide")
st.title("🧬 PCOS Pre-Ranked GSEA Dashboard")
st.markdown("Transcriptomic analysis of Polycystic Ovary Syndrome (PCOS) evaluating pathway-level shifts via Gene Ontology (GO).")

# Cache the data loading and GSEA so it doesn't re-run on every click
@st.cache_data
def run_gsea_pipeline():
    # Load data
    results_df = pd.read_csv("pcos_deseq2_results.csv", index_col=0)
    annot_df = pd.read_csv("Human.GRCh38.p13.annot.tsv", sep='\t', dtype=str)
    
    # Merge annotations
    results_df = results_df.reset_index()
    results_df['GeneID'] = results_df['GeneID'].astype(str)
    results_df = results_df.merge(annot_df[['GeneID', 'Symbol']], on='GeneID', how='left')
    results_df['Gene_Name'] = results_df['Symbol'].fillna(results_df['GeneID'])
    
    # Clean and rank
    clean_df = results_df.dropna(subset=['Gene_Name', 'log2FoldChange'])
    rnk_df = clean_df[['Gene_Name', 'log2FoldChange']].sort_values(by='log2FoldChange', ascending=False)
    
    # Run GSEA
    pre_res = gp.prerank(
        rnk=rnk_df,
        gene_sets='GO_Biological_Process_2021',
        threads=4,
        min_size=5,
        max_size=1000,
        permutation_num=1000,
        outdir=None,
        seed=42
    )
    return pre_res.res2d

with st.spinner("Running 1000 permutations for GSEA... This takes ~15 seconds on the first load."):
    gsea_results = run_gsea_pipeline()

# Filter for significance (FDR < 0.05) and grab top 15 by absolute NES
sig_results = gsea_results[gsea_results['FDR q-val'] < 0.05].copy()
sig_results['Abs_NES'] = sig_results['NES'].abs()
top_pathways = sig_results.sort_values(by='Abs_NES', ascending=False).head(15)

# Visualizations
st.subheader("📊 Top Enriched Biological Processes (FDR < 0.05)")

# Create a Plotly Bar Chart
fig = px.bar(
    top_pathways,
    x="NES",
    y="Term",
    orientation='h',
    color="FDR q-val",
    color_continuous_scale="Viridis",
    title="Normalized Enrichment Score (NES) by Pathway"
)
fig.update_layout(yaxis={'categoryorder':'total ascending'}, height=500)
st.plotly_chart(fig, use_container_width=True)

# Data Table
st.subheader("📋 Full GSEA Results Inspector")
st.dataframe(gsea_results[['Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes']], use_container_width=True)
