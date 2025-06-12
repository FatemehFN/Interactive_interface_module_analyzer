import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import scanpy as sc




# Optional: Tweak font size to avoid warnings from scanpy
matplotlib.rcParams.update({'font.size': 10})
from backend import (
    build_knn_graph,
    run_leiden,
    compute_eigengenes,
    correlate_eigengenes,
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
                corr = correlate_eigengenes(eigengenes_df, pheno_df)
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

                sns.heatmap(
                    corr.T,  # keep labels by using DataFrame
                    annot=True,
                    fmt=".2f",
                    cmap="vlag",
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={"label": "Pearson r"},
                    annot_kws={"size": font_scale * 10},
                    ax=ax
                )


                plt.xlabel("Modules")
                plt.ylabel("Phenotypes")
                st.pyplot(fig)



                st.subheader("🧠 Gene Set Enrichment by Cell Type")
                enrichment_df, enrichment_fig = perform_cell_type_enrichment_and_heatmap(expr_df, clusters.tolist())
                st.dataframe(enrichment_df)
                st.pyplot(enrichment_fig)


            except Exception as e:
                st.error(f"❌ An error occurred during analysis: {e}")



# if expr_df is not None:

#     if st.button("Analyze k"):
#         with st.spinner("Running k analysis..."):
#             k_range = list(range(2, 10))
#             best_k_df = analyze_best_k(expr_df, k_range, resolution)
#             st.dataframe(best_k_df)
#             fig2, ax2 = plt.subplots(figsize=(10, 5))
#             ax2.plot(best_k_df['k'], best_k_df['n_clusters'], marker='o')
#             ax2.set_xlabel('k (neighbors)')
#             ax2.set_ylabel('Number of Leiden clusters')
#             ax2.set_title('Number of clusters vs k (res=1.0)')
#             st.pyplot(fig2)