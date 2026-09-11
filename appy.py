import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx
import requests

st.set_page_config(page_title="PCOS Transcriptomics", layout="wide")
st.title("🧬 PCOS Transcriptomic Profiling Dashboard")
st.markdown("Transcriptomic analysis of Polycystic Ovary Syndrome (PCOS) evaluating differential expression, functional pathways, and master regulators using strictly validated mRNA RNA-seq samples.")

# --- DATA LOADING ---
@st.cache_data
def load_data():
    # 1. Load results and raw counts
    results_df = pd.read_csv("pcos_deseq2_results.csv", index_col=0).reset_index()
    annot_df = pd.read_csv("Human.GRCh38.p13.annot.tsv", sep='\t', dtype=str)
    gsea_df = pd.read_csv("gsea_results_precalculated.csv")
    tf_df = pd.read_csv("tf_screening_results.csv")
    raw_counts = pd.read_csv("GSE168404_mRNA_C-vs-P.all.txt", index_col=0, sep='\t', encoding='utf-16')
    
    # 2. Merge gene names safely
    results_df['GeneID'] = results_df['GeneID'].astype(str)
    results_df = results_df.merge(annot_df[['GeneID', 'Symbol']], on='GeneID', how='left')
    results_df['Gene_Name'] = results_df['Symbol'].fillna(results_df['GeneID'])
    
    # 3. Prepare DESeq2 DataFrame
    de_df = results_df.dropna(subset=['padj', 'log2FoldChange', 'baseMean']).copy()
    de_df['-log10(padj)'] = -np.log10(de_df['padj'] + 1e-300)

    def categorize(row):
        if row['padj'] < 0.05 and row['log2FoldChange'] > 1: return 'Upregulated'
        if row['padj'] < 0.05 and row['log2FoldChange'] < -1: return 'Downregulated'
        return 'Not Significant'

    de_df['Significance'] = de_df.apply(categorize, axis=1)
    
    return de_df, gsea_df, tf_df, raw_counts

de_df, gsea_df, tf_df, raw_counts = load_data()
color_map = {'Upregulated': '#EF553B', 'Downregulated': '#636EFA', 'Not Significant': '#E5E5E5'}

# --- STRING API NETWORK FETCH ---
@st.cache_data
def fetch_ppi_network(gsea_dataframe):
    target_terms = gsea_dataframe.sort_values(by='NES', key=abs, ascending=False).head(2)['Term'].tolist()
    genes = []
    for term in target_terms:
        gene_string = gsea_dataframe[gsea_dataframe['Term'] == term]['Lead_genes'].iloc[0]
        if pd.notna(gene_string):
            genes.extend(str(gene_string).split(';')[:15])
    
    genes = list(set(genes))
    if not genes:
        return None
        
    string_api_url = "https://string-db.org/api/json/network"
    params = {"identifiers": "%0d".join(genes), "species": 9606, "caller_identity": "pcos_dashboard"}
    
    try:
        response = requests.post(string_api_url, data=params, timeout=15)
        response.raise_for_status()
        return response.json()
    except:
        return None

# --- DASHBOARD TABS ---
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Differential Expression", 
    "🧬 Pathways (GSEA)", 
    "⚙️ Master Regulators", 
    "🕸️ PPI Network",
    "📋 Data Inspector"
])

# --- TAB 1: DESeq2 GENE LEVEL ---
with tab1:
    st.subheader("🌋 Volcano Plot")
    fig_volcano = px.scatter(
        de_df, x='log2FoldChange', y='-log10(padj)', color='Significance',
        color_discrete_map=color_map, hover_name='Gene_Name', height=500
    )
    fig_volcano.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black")
    fig_volcano.add_vline(x=1, line_dash="dash", line_color="black")
    fig_volcano.add_vline(x=-1, line_dash="dash", line_color="black")
    st.plotly_chart(fig_volcano, use_container_width=True)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("📈 MA Plot")
        fig_ma = px.scatter(
            de_df, x='baseMean', y='log2FoldChange', color='Significance', log_x=True,
            color_discrete_map=color_map, hover_name='Gene_Name'
        )
        fig_ma.add_hline(y=0, line_dash="dash", line_color="black")
        st.plotly_chart(fig_ma, use_container_width=True)

    with col2:
        st.subheader("🔥 Top 40 Genes Heatmap")
        valid_samples = ['C1', 'C2', 'C3', 'C4', 'C5', 'P1', 'P2', 'P3', 'P4', 'P5']
        sig_genes = de_df.sort_values('padj').head(40)
        
        # Safe Pandas Z-Score calculation
        heat_df = raw_counts.loc[sig_genes['GeneID'], valid_samples].copy().astype(float)
        heat_df.index = sig_genes['Gene_Name'].values
        heat_mean = heat_df.mean(axis=1)
        heat_std = heat_df.std(axis=1).replace(0, 1e-9)
        heat_z = heat_df.sub(heat_mean, axis=0).div(heat_std, axis=0)

        fig_heat = px.imshow(
            heat_z, x=valid_samples, y=heat_z.index, color_continuous_scale="RdBu_r", aspect="auto"
        )
        fig_heat.update_xaxes(side="top")
        st.plotly_chart(fig_heat, use_container_width=True)

# --- TAB 2: GSEA PATHWAY LEVEL ---
with tab2:
    st.subheader("📊 Top Enriched Biological Processes")
    sig_results = gsea_df[gsea_df['FDR q-val'] < 0.05].copy()
    sig_results['Abs_NES'] = sig_results['NES'].abs()
    top_pathways = sig_results.sort_values(by='Abs_NES', ascending=False).head(15)

    fig_bar = px.bar(
        top_pathways, x="NES", y="Term", orientation='h', color="FDR q-val", color_continuous_scale="Viridis"
    )
    fig_bar.update_layout(yaxis={'categoryorder':'total ascending'}, height=450)
    st.plotly_chart(fig_bar, use_container_width=True)

    st.subheader("🫧 Pathway Bubble Plot")
    top_pathways['Gene_Count'] = top_pathways['Lead_genes'].apply(lambda x: len(str(x).split(';')))
    fig_bubble = px.scatter(
        top_pathways, x='NES', y='Term', size='Gene_Count', color='FDR q-val', color_continuous_scale='Viridis_r',
        hover_data=['Gene_Count', 'FDR q-val']
    )
    fig_bubble.update_layout(yaxis={'categoryorder':'total ascending'}, height=450)
    st.plotly_chart(fig_bubble, use_container_width=True)

# --- TAB 3: MASTER REGULATORS (TFs) ---
with tab3:
    st.subheader("⚙️ Upstream Transcription Factor Screening (ChIP-X)")
    sig_tfs = tf_df[tf_df['FDR q-val'] < 0.05].copy()
    sig_tfs['Abs_NES'] = sig_tfs['NES'].abs()
    top_tfs = sig_tfs.sort_values(by='Abs_NES', ascending=False).head(15)

    fig_tf = px.bar(
        top_tfs, x="NES", y="Term", orientation='h', color="FDR q-val", color_continuous_scale="Cividis"
    )
    fig_tf.update_layout(yaxis={'categoryorder':'total ascending'}, height=500)
    st.plotly_chart(fig_tf, use_container_width=True)

# --- TAB 4: PPI NETWORK ---
with tab4:
    st.subheader("🕸️ Protein-Protein Interaction (PPI) Network")
    st.markdown("Network mapping of the leading-edge core drivers from the top two dysregulated pathways.")
    
    interactions = fetch_ppi_network(gsea_df)
    
    if interactions:
        G = nx.Graph()
        for edge in interactions:
            if edge['score'] > 0.4:
                G.add_edge(edge['preferredName_A'], edge['preferredName_B'], weight=edge['score'])

        if len(G.nodes) > 0:
            pos = nx.spring_layout(G, seed=42, k=0.5)
            edge_x, edge_y = [], []
            for edge in G.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

            edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=1, color='#888'), hoverinfo='none', mode='lines')

            node_x = [pos[node][0] for node in G.nodes()]
            node_y = [pos[node][1] for node in G.nodes()]
            node_text = [node for node in G.nodes()]
            node_adjacencies = [len(list(G.neighbors(node))) for node in G.nodes()]

            node_trace = go.Scatter(
                x=node_x, y=node_y, mode='markers+text', text=node_text, textposition="top center", hoverinfo='text',
                marker=dict(
                    showscale=True, colorscale='Viridis', size=15, color=node_adjacencies,
                    colorbar=dict(title="Connections", thickness=15, xanchor='left'), line_width=2
                )
            )

            fig_net = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(showlegend=False, hovermode='closest', margin=dict(b=0,l=0,r=0,t=0),
                                 xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                 yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), height=600)
            )
            st.plotly_chart(fig_net, use_container_width=True)
        else:
            st.warning("No high-confidence connections found among the top driver genes.")
    else:
        st.warning("Could not fetch network data. Please check your internet connection or try again later.")

# --- TAB 5: DATA INSPECTOR ---
with tab5:
    st.subheader("📋 DESeq2 Differential Expression Table")
    st.dataframe(de_df[['GeneID', 'Gene_Name', 'baseMean', 'log2FoldChange', 'padj', 'Significance']], use_container_width=True)
    
    st.markdown("---")
    st.subheader("📋 GSEA Pathway Results")
    st.dataframe(gsea_df[['Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes']], use_container_width=True)
