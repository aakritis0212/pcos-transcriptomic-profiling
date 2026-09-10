import streamlit as st
import pandas as pd
import numpy as np
import gseapy as gp
import plotly.express as px

st.set_page_config(page_title="PCOS Transcriptomics", layout="wide")
st.title("🧬 PCOS Pre-Ranked GSEA Dashboard")
st.markdown("Transcriptomic analysis of Polycystic Ovary Syndrome (PCOS) evaluating pathway-level shifts via Gene Ontology (GO).")

# 1. Load Data (Cached so it only runs once)
@st.cache_data
def load_data():
    df = pd.read_csv("pcos_deseq2_results.csv", index_col=0)
    annot = pd.read_csv("Human.GRCh38.p13.annot.tsv", sep='\t', dtype=str)
    
    # Merge annotations
    df = df.reset_index()
    df['GeneID'] = df['GeneID'].astype(str)
    df = df.merge(annot[['GeneID', 'Symbol']], on='GeneID', how='left')
    df['Gene_Name'] = df['Symbol'].fillna(df['GeneID'])
    return df

# 2. Run GSEA (Cached so it doesn't re-run on every click)
@st.cache_data
def run_gsea_pipeline(df):
    clean_df = df.dropna(subset=['Gene_Name', 'log2FoldChange'])
    rnk_df = clean_df[['Gene_Name', 'log2FoldChange']].sort_values(by='log2FoldChange', ascending=False)
    
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

# --- EXECUTE PIPELINE ---
results_df = load_data()

with st.spinner("Running 1000 permutations for GSEA... This takes ~15 seconds on the first load."):
    gsea_results = run_gsea_pipeline(results_df)

# --- VISUALIZATION 1: GSEA BAR CHART ---
st.subheader("📊 Top Enriched Biological Processes (FDR < 0.05)")

sig_results = gsea_results[gsea_results['FDR q-val'] < 0.05].copy()
sig_results['Abs_NES'] = sig_results['NES'].abs()
top_pathways = sig_results.sort_values(by='Abs_NES', ascending=False).head(15)

fig_bar = px.bar(
    top_pathways,
    x="NES",
    y="Term",
    orientation='h',
    color="FDR q-val",
    color_continuous_scale="Viridis",
    title="Normalized Enrichment Score (NES) by Pathway"
)
fig_bar.update_layout(yaxis={'categoryorder':'total ascending'}, height=500)
st.plotly_chart(fig_bar, use_container_width=True)

# --- VISUALIZATION 2: VOLCANO PLOT ---
st.markdown("---")
st.subheader("🌋 Differential Expression: Volcano Plot")

# Clean the data for the plot
volcano_df = results_df.dropna(subset=['padj', 'log2FoldChange']).copy()
volcano_df['-log10(padj)'] = -np.log10(volcano_df['padj'] + 1e-300)

def categorize_significance(row):
    if row['padj'] < 0.05 and row['log2FoldChange'] > 1:
        return 'Upregulated'
    elif row['padj'] < 0.05 and row['log2FoldChange'] < -1:
        return 'Downregulated'
    else:
        return 'Not Significant'

volcano_df['Significance'] = volcano_df.apply(categorize_significance, axis=1)

fig_volcano = px.scatter(
    volcano_df,
    x='log2FoldChange',
    y='-log10(padj)',
    color='Significance',
    color_discrete_map={
        'Upregulated': '#EF553B',
        'Downregulated': '#636EFA',
        'Not Significant': '#E5E5E5'
    },
    hover_name=volcano_df['Gene_Name'], 
    labels={
        'log2FoldChange': 'Log2 Fold Change',
        '-log10(padj)': '-Log10(Adjusted P-value)'
    }
)

# Add threshold lines
fig_volcano.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black")
fig_volcano.add_vline(x=1, line_dash="dash", line_color="black")
fig_volcano.add_vline(x=-1, line_dash="dash", line_color="black")

st.plotly_chart(fig_volcano, use_container_width=True)

# --- DATA TABLE ---
st.markdown("---")
st.subheader("📋 Full GSEA Results Inspector")
st.dataframe(gsea_results[['Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes']], use_container_width=True)
