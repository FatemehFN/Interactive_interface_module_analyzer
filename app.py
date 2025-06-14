import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import scanpy as sc
import numpy as np


# Optional: Tweak font size to avoid warnings from scanpy
matplotlib.rcParams.update({'font.size': 10})
from backend import (
    build_knn_graph,
    run_leiden,
    compute_eigengenes,
    correlate_eigengenes,
    correlate_eigengenes_Biweight_midcorrelation,
    perform_cell_type_enrichment_and_heatmap,
    run_infomap
)

st.set_page_config(page_title="Network Module Analyzer", layout="wide")
st.title("🔬 Network Module Analyzer")

st.markdown("""
Upload your expression matrix (genes × samples) and phenotype data (samples × traits),  
or use our example dataset to explore the analysis pipeline.
""")

# --- Upload Option ---
use_own_data = st.radio(
    "📂 Do you want to upload your own data?",
    ("No, use example data", "Yes, I want to upload my own"),
    index=0
)

# --- Load Data ---
expr_df, pheno_df = None, None

if use_own_data == "Yes, I want to upload my own":
    expr_file = st.file_uploader("🧬 Upload Expression Data (CSV)", type="csv")
    pheno_file = st.file_uploader("📊 Upload Phenotype Data (CSV)", type="csv")

    if expr_file and pheno_file:
        expr_df = pd.read_csv(expr_file, index_col=0)
        pheno_df = pd.read_csv(pheno_file, index_col=0)
else:
    try:
        expr_df = pd.read_csv("example_data/expression.csv", index_col=0)
        pheno_df = pd.read_csv("example_data/phenotype.csv", index_col=0)
        st.success("✅ Example data loaded.")
    except FileNotFoundError:
        st.error("❌ Example data not found. Please make sure `example_data/expression.csv` and `example_data/phenotype.csv` exist.")
        st.stop()

# --- Parameter Selection ---
k = st.slider("🔧 Select `k` for KNN graph", min_value=2, max_value=50, value=10)
# --- Clustering method selection ---
clustering_method = st.radio("🤖 Choose clustering algorithm:", options=["Leiden", "Infomap"], index=0)

if clustering_method == "Leiden":
    resolution = st.slider("🌀 Select resolution for Leiden clustering", min_value=0.1, max_value=2.0, value=1.0, step=0.1)
else:
    markov_time = st.slider("⏱️ Select Markov time for Infomap", min_value=0.1, max_value=5.0, value=1.0, step=0.1)

if expr_df is not None and pheno_df is not None:
    st.write("Expression Data:", expr_df.shape)
    st.write("Phenotype Data:", pheno_df.shape)

    # --- Correlation method selection ---
    correlation_method = st.radio("📊 Choose correlation method:", options=["Pearson", "Biweight midcorrelation"], index=0)

    if st.button("🚀 Run Analysis"):
        with st.spinner("Processing..."):
            try:
                # Step 1: KNN Graph and Clustering
                adata = build_knn_graph(expr_df, k)

                if clustering_method == "Leiden":
                    clusters = run_leiden(adata, resolution)
                else:
                    clusters = run_infomap(adata, markov_time)

                st.subheader("📂 Clusters")
                st.dataframe(clusters.value_counts().sort_index())

                # --- UMAP + t-SNE side-by-side ---
                st.subheader("🧭 UMAP and 🌀 t-SNE of Samples")

                fig, axes = plt.subplots(1, 2, figsize=(14, 6))

                
                if clustering_method == "Leiden":
                    # UMAP
                    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='UMAP')
                    # t-SNE
                    sc.pl.tsne(adata, color='leiden', ax=axes[1], show=False, title='t-SNE')
                
                elif clustering_method == "Infomap":
                    # UMAP
                    sc.pl.umap(adata, color='infomap', ax=axes[0], show=False, title='UMAP')
                    # t-SNE
                    sc.pl.tsne(adata, color='infomap', ax=axes[1], show=False, title='t-SNE')

                plt.tight_layout()
                st.pyplot(fig)

                # Step 2: Eigengenes
                st.subheader("🧬 Module Eigengenes")
                st.write("Expression matrix shape:", expr_df.shape)
                eigengenes_df = compute_eigengenes(expr_df, clusters.tolist())
                st.dataframe(eigengenes_df)

                # Step 3: Correlation
                st.subheader("📈 Correlation with Phenotypes")
                if correlation_method == "Pearson":
                    corr, pval = correlate_eigengenes(eigengenes_df, pheno_df)
                elif correlation_method == "Biweight midcorrelation":
                     corr, pval = correlate_eigengenes_Biweight_midcorrelation(eigengenes_df, pheno_df)
                st.dataframe(corr)

                corr=corr.astype(float)
                st.subheader("📊 Heatmap of Correlations (Phenotypes x Modules)")
                # Transpose if needed: rows = phenotypes, columns = modules
                # Dynamically scale figure height
                fig_height = max(4, len(corr.columns))  # based on number of phenotypes
                fig, ax = plt.subplots(figsize=(10, fig_height))

                # Dynamically adjust font size
                n_modules = corr.shape[0]
                n_phenotypes = corr.shape[1]
                font_scale = min(1.2, max(0.4, 30 / (n_modules * n_phenotypes)))

                # Calculate symmetric limits around 0
                vmax = max(abs(corr.min().min()), abs(corr.max().max()))
                vmin = -vmax

                # Draw the heatmap with correlation values
                sns.heatmap(
                    corr.T,  # keep labels by using DataFrame
                    annot=True, # Display numerical values
                    fmt=".2f",
                    cmap="vlag",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"label": f"{correlation_method} r"}, # Use selected method in label
                    annot_kws={"size": font_scale * 14}, # Increased font size for correlation numbers (from 12 to 14)
                    ax=ax
                )

                # Add stars for p < 0.05 on top of correlation numbers
                for i in range(pval.T.shape[0]): # Iterate over rows of pval.T (phenotypes)
                    for j in range(pval.T.shape[1]): # Iterate over columns of pval.T (modules)
                        if pval.T.iloc[i, j] < 0.05:
                            # Add a text annotation for the star
                            # Adjusted y-coordinate to move the star further up (0.5 - 0.2 -> 0.5 - 0.3)
                            ax.text(j + 0.5, i + 0.5 - 0.25, '*', ha='center', va='center',
                                    color='black', fontsize=font_scale * 20, weight='bold') # Adjusted va to 'center' and moved using y-offset


                plt.xlabel("Modules")
                plt.ylabel("Phenotypes")
                st.pyplot(fig)


                st.subheader("🧠 Gene Set Enrichment by Cell Type")
                enrichment_df, enrichment_fig = perform_cell_type_enrichment_and_heatmap(expr_df, clusters.tolist())
                st.dataframe(enrichment_df)
                st.pyplot(enrichment_fig)


            except Exception as e:
                st.error(f"❌ An error occurred during analysis: {e}")
