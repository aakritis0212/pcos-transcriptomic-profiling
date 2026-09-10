import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

st.set_page_config(page_title="PCOS Transcriptomics", layout="wide")
st.title("🧬 PCOS Transcriptomic Profiling & Enrichment Dashboard")
st.markdown("Transcriptomic analysis of Polycystic Ovary Syndrome (PCOS) evaluating differential expression and Gene Ontology (GO) pathway shifts.")

# Load Cached Data
@st.cache_data
def load_data():
    results_df = pd.read_csv("pcos_deseq2_results.csv", index_col=0)
    annot_df = pd.read_csv("Human.GRCh38.p13.annot.tsv", sep='\t', dtype=str)
    
    results_df = results_df.reset_index()
    results_df['GeneID'] = results_df['GeneID'].astype(str)
    results_df = results_df.merge(annot_df[['GeneID', 'Symbol']], on='GeneID', how='left')
    results_df['Gene_Name'] = results_df['Symbol'].fillna(results_df['GeneID'])
    
    gsea_results = pd.read_csv("gsea_results_precalculated.csv")
    return results_df, gsea_results

results_df, gsea_results = load_data()

# Prepare clean DESeq2 data
de_df = results_df.dropna(subset=['padj', 'log2FoldChange', 'baseMean']).copy()
de_df['-log10(padj)'] = -np.log10(de_df['padj'] + 1e-300)

def categorize_significance(row):
    if row['padj'] < 0.05 and row['log2FoldChange'] > 1:
        return 'Upregulated'
    elif row['padj'] < 0.05 and row['log2FoldChange'] < -1:
        return 'Downregulated'
    else:
        return 'Not Significant'

de_df['Significance'] = de_df.apply(categorize_significance, axis=1)

# Tabbed Interface
tab1, tab2, tab3 = st.tabs(["📊 DESeq2 Gene Level", "🧬 GSEA Pathway Level", "📋 Data Inspector"])

# --- TAB 1: DESeq2 GENE LEVEL ---
with tab1:
    st.subheader("🌋 Volcano Plot")
    fig_volcano = px.scatter(
        de_df,
        x='log2FoldChange',
        y='-log10(padj)',
        color='Significance',
        render_mode='webgl',
        color_discrete_map={
            'Upregulated': '#EF553B',
            'Downregulated': '#636EFA',
            'Not Significant': '#E5E5E5'
        },
        hover_name='Gene_Name',
        labels={'log2FoldChange': 'Log2 Fold Change', '-log10(padj)': '-Log10(Adjusted P-value)'}
    )
    fig_volcano.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black")
    fig_volcano.add_vline(x=1, line_dash="dash", line_color="black")
    fig_volcano.add_vline(x=-1, line_dash="dash", line_color="black")
    st.plotly_chart(fig_volcano, use_container_width=True)

    st.markdown("---")
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📈 MA Plot")
        fig_ma = px.scatter(
            de_df,
            x='baseMean',
            y='log2FoldChange',
            color='Significance',
            log_x=True,
            render_mode='webgl',
            color_discrete_map={
                'Upregulated': '#EF553B',
                'Downregulated': '#636EFA',
                'Not Significant': '#E5E5E5'
            },
            hover_name='Gene_Name',
            labels={'baseMean': 'Mean Expression (baseMean, Log Scale)', 'log2FoldChange': 'Log2 Fold Change'}
        )
        fig_ma.add_hline(y=0, line_dash="dash", line_color="black")
        st.plotly_chart(fig_ma, use_container_width=True)

    with col2:
        st.subheader("🏆 Top 20 Differentially Expressed Genes")
        sig_genes = de_df[de_df['padj'] < 0.05].copy()
        top_20 = sig_genes.sort_values(by='padj', ascending=True).head(20).sort_values(by='log2FoldChange', ascending=True)
        
        fig_top = px.bar(
            top_20,
            x='log2FoldChange',
            y='Gene_Name',
            orientation='h',
            color='log2FoldChange',
            color_continuous_scale='Bluered',
            labels={'log2FoldChange': 'Log2 Fold Change', 'Gene_Name': 'Gene'}
        )
        fig_top.update_layout(yaxis={'categoryorder':'total ascending'})
        st.plotly_chart(fig_top, use_container_width=True)

# --- TAB 2: GSEA PATHWAY LEVEL ---
with tab2:
    sig_results = gsea_results[gsea_results['FDR q-val'] < 0.05].copy()
    sig_results['Abs_NES'] = sig_results['NES'].abs()
    top_pathways = sig_results.sort_values(by='Abs_NES', ascending=False).head(15)

    st.subheader("📊 Top Enriched Biological Processes (NES)")
    fig_bar = px.bar(
        top_pathways,
        x="NES",
        y="Term",
        orientation='h',
        color="FDR q-val",
        color_continuous_scale="Viridis",
        title="Normalized Enrichment Score (NES) by Pathway"
    )
    fig_bar.update_layout(yaxis={'categoryorder':'total ascending'}, height=450)
    st.plotly_chart(fig_bar, use_container_width=True)

    st.markdown("---")
    st.subheader("🫧 GSEA Pathway Bubble Plot")
    sig_results['Gene_Count'] = sig_results['Lead_genes'].apply(lambda x: len(str(x).split(';')))
    top_bubbles = sig_results.sort_values(by='Abs_NES', ascending=False).head(15)

    fig_bubble = px.scatter(
        top_bubbles,
        x='NES',
        y='Term',
        size='Gene_Count',
        color='FDR q-val',
        color_continuous_scale='Viridis_r',
        hover_data=['Gene_Count', 'FDR q-val'],
        labels={'NES': 'Normalized Enrichment Score', 'Term': 'Biological Process'}
    )
    fig_bubble.update_layout(yaxis={'categoryorder':'total ascending'}, height=450)
    st.plotly_chart(fig_bubble, use_container_width=True)

# --- TAB 3: DATA INSPECTOR ---
with tab3:
    st.subheader("📋 GSEA Pathway Results")
    st.dataframe(gsea_results[['Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes']], use_container_width=True)
    
    st.markdown("---")
    st.subheader("📋 DESeq2 Differential Expression Table")
    st.dataframe(de_df[['GeneID', 'Gene_Name', 'baseMean', 'log2FoldChange', 'pvalue', 'padj', 'Significance']], use_container_width=True)
