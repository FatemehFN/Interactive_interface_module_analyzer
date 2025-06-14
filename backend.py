import pandas as pd
import scanpy as sc
import numpy as np
from sklearn.decomposition import PCA
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import fisher_exact
import seaborn as sns
import matplotlib.pyplot as plt
import io
import base64
import json
from scipy.stats import fisher_exact
from infomap import Infomap
import networkx as nx
import pingouin as pg




def build_knn_graph(expr_matrix, k):
    adata = sc.AnnData(expr_matrix)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=k)
    
    # Compute UMAP and t-SNE for visualization
    sc.tl.umap(adata)
    sc.tl.tsne(adata, n_pcs=20)  # You can adjust this
    return adata

def run_leiden(adata, resolution):
    sc.tl.leiden(adata, resolution=resolution)
    return adata.obs['leiden']


def compute_eigengenes(expr_matrix, cluster_labels):
    """
    Compute module eigengenes using the first principal component of each module,
    and ensure consistent sign by aligning eigengene with first gene in the module.

    Parameters:
        expr_matrix (pd.DataFrame): genes × samples
        cluster_labels (list or pd.Series): cluster labels for each gene (length = number of genes)

    Returns:
        pd.DataFrame: samples × modules (Module_0, Module_1, ...)
    """
    if len(cluster_labels) != expr_matrix.shape[0]:
        raise ValueError("Length of cluster_labels must match number of genes (rows in expr_matrix).")

    modules = {}

    try:
        unique_labels = sorted(set(cluster_labels), key=lambda x: int(x))
    except ValueError:
        unique_labels = sorted(set(cluster_labels))

    for label in unique_labels:
        idx = [i for i, val in enumerate(cluster_labels) if val == label]
        sub_expr = expr_matrix.iloc[idx, :]  # genes × samples

        pca = PCA(n_components=1)
        eigengene = pca.fit_transform(sub_expr.T).flatten()  # 1 value per sample

        # Align sign with the first gene in the module
        first_gene_expr = sub_expr.iloc[0, :].values
        sign = np.sign(np.corrcoef(eigengene, first_gene_expr)[0, 1])
        eigengene *= sign

        modules[f"Module_{label}"] = eigengene

    return pd.DataFrame(modules, index=expr_matrix.columns)



from scipy.stats import pearsonr
import pandas as pd
import numpy as np

def correlate_eigengenes(eigengenes, phenotypes):
    """
    Compute Pearson correlations and p-values between module eigengenes and phenotype traits.

    Parameters:
        eigengenes (pd.DataFrame): samples × modules
        phenotypes (pd.DataFrame): samples × traits

    Returns:
        (pd.DataFrame, pd.DataFrame): correlation matrix and p-value matrix (modules × traits)
    """
    # Ensure indices are strings for alignment
    eigengenes.index = eigengenes.index.astype(str)
    phenotypes.index = phenotypes.index.astype(str)

    # Align samples
    common_samples = eigengenes.index.intersection(phenotypes.index)
    eigengenes = eigengenes.loc[common_samples]
    phenotypes = phenotypes.loc[common_samples]

    # Ensure numeric types
    eigengenes = eigengenes.apply(pd.to_numeric, errors='coerce')
    phenotypes = phenotypes.apply(pd.to_numeric, errors='coerce')

    # Drop samples with NaNs in either
    aligned = eigengenes.join(phenotypes, how='inner')
    eigengenes = aligned[eigengenes.columns]
    phenotypes = aligned[phenotypes.columns]

    # Initialize result DataFrames
    corr_df = pd.DataFrame(index=eigengenes.columns, columns=phenotypes.columns, dtype=float)
    pval_df = pd.DataFrame(index=eigengenes.columns, columns=phenotypes.columns, dtype=float)

    # Compute correlation and p-values
    for m in eigengenes.columns:
        for ph in phenotypes.columns:
            # Drop NaNs for this pair
            df = pd.DataFrame({'x': eigengenes[m], 'y': phenotypes[ph]})
            df = df.replace([np.inf, -np.inf], np.nan).dropna()

            if len(df) >= 2:  # Need at least two samples for Pearson
                r, p = pearsonr(df['x'], df['y'])
                corr_df.loc[m, ph] = r
                pval_df.loc[m, ph] = p
            else:
                corr_df.loc[m, ph] = np.nan
                pval_df.loc[m, ph] = np.nan

    return corr_df, pval_df




def correlate_eigengenes_Biweight_midcorrelation(eigengenes, phenotypes):
    """
    Compute Pearson correlations between module eigengenes and phenotype traits.

    Parameters:
        eigengenes (pd.DataFrame): samples × modules
        phenotypes (pd.DataFrame): samples × traits

    Returns:
        pd.DataFrame: modules × traits correlation matrix
    """

    # Ensure both indices are strings for proper alignment
    eigengenes.index = eigengenes.index.astype(str)
    phenotypes.index = phenotypes.index.astype(str)

    # Align samples (rows)
    common_samples = eigengenes.index.intersection(phenotypes.index)
    #print(f"Found {len(common_samples)} common samples.")

    eigengenes = eigengenes.loc[common_samples]
    phenotypes = phenotypes.loc[common_samples]

    # Convert to numeric in case of non-numeric columns
    eigengenes = eigengenes.apply(pd.to_numeric, errors='coerce')
    phenotypes = phenotypes.apply(pd.to_numeric, errors='coerce')

    # Drop samples with NaNs in either dataframe
    aligned = eigengenes.join(phenotypes, how='inner')
    eigengenes = aligned[eigengenes.columns]
    phenotypes = aligned[phenotypes.columns]

    # Initialize correlation DataFrame
    corr_df = pd.DataFrame(index=eigengenes.columns, columns=phenotypes.columns, dtype=float)
    pvalue_df = pd.DataFrame(index=eigengenes.columns, columns=phenotypes.columns, dtype=float)

    # Compute correlations
    for m in eigengenes.columns:
        for ph in phenotypes.columns:

            result = pg.corr(eigengenes[m], phenotypes[ph], method='bicor')
            corr_df.loc[m, ph] = result['r'].values[0]
            pvalue_df.loc[m, ph] = result['p-val'].values[0]

    return corr_df, pvalue_df






# def analyze_best_k(expr_matrix, k_values, resolution=1.0):
#     results = []
#     for k in k_values:
#         adata = build_knn_graph(expr_matrix, k)
#         sc.tl.leiden(adata, resolution=resolution)
#         clusters = adata.obs['leiden']
        
#         n_clusters = clusters.nunique()
#         # Optional: compute modularity (scanpy has it in adata.uns['modularity'] if computed)
#         # or silhouette score - but silhouette needs embedding or distance matrix
        
#         results.append({'k': k, 'n_clusters': n_clusters})
#     return pd.DataFrame(results)



def perform_cell_type_enrichment_and_heatmap(expr_df, clusters, threshold=1e-3):

    # Define marker genes
    brain_marker_genes = {
        "Neurons": ["MAP2", "NEFL", "SYT1", "GRIN2D", "SNAP23", "DRD3", "SLC6A13", "GNAO1", "MEF2C", "CNR1", "TUBB4A"],
        "Astrocytes": ["GFAP", "S100B", "AQP4", "ALDH1A1", "SLC1A3", "GLUL", "CD44", "GJA1", "VIM"],
        "Oligodendrocytes": ["MBP", "PLP1", "MOG", "OLIG1", "OLIG2", "SOX10", "NKX6-2", "CNP"],
        "Microglia": ["AIF1", "CD68", "TMEM119", "ITGAM", "TREM2", "P2RY12", "CX3CR1", "SPI1", "CD14"],
        "Endothelial Cells": ["CLDN5", "PECAM1", "FLT1", "ENG", "CD34", "VWF", "ICAM1"],
        "Inhibitory Neurons": ["GAD1", "SLC32A1", "SST", "PVALB"],
        "Excitatory Neurons": ["SLC17A7", "GRIN2A", "SYT1"],
        "Dopaminergic Neurons": ["TH", "DRD1", "DRD2", "DBH"],
        "Cholinergic Neurons": ["CHAT", "AChE"]
    }

    # Map genes to clusters
    gene_clusters = pd.Series(clusters, index=expr_df.index)
    unique_modules = sorted(gene_clusters.unique(), key=int)

    enrichment_matrix = pd.DataFrame(index=unique_modules, columns=brain_marker_genes.keys())
    pval_matrix = pd.DataFrame(index=unique_modules, columns=brain_marker_genes.keys())

    universe_genes = set(expr_df.index)

    for module in enrichment_matrix.index:
        module_genes = set(gene_clusters[gene_clusters == module].index)

        for cell_type, markers in brain_marker_genes.items():
            marker_set = set(markers) & universe_genes  # markers in dataset

            # Fisher's exact test
            a = len(module_genes & marker_set)
            b = len(module_genes) - a
            c = len(marker_set) - a
            d = len(universe_genes) - a - b - c

            _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
            enrichment_score = -np.log10(p_value + 1e-10)

            enrichment_matrix.loc[module, cell_type] = enrichment_score
            pval_matrix.loc[module, cell_type] = p_value

    # Transpose for heatmap: cell types × modules
    enrichment_matrix = enrichment_matrix.astype(float).T
    pval_matrix = pval_matrix.astype(float).T

    # Generate star annotations (above number if p < 0.05)
    annot_matrix = pd.DataFrame(index=enrichment_matrix.index, columns=enrichment_matrix.columns)
    for row in enrichment_matrix.index:
        for col in enrichment_matrix.columns:
            score = enrichment_matrix.loc[row, col]
            pval = pval_matrix.loc[row, col]
            if score < threshold:
                annot_matrix.loc[row, col] = ""
            else:
                star = "*" if pval < 0.05 else ""
                annot_matrix.loc[row, col] = f"{star}\n{score:.2f}"

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        enrichment_matrix,
        annot=annot_matrix,
        fmt="",
        cmap="YlGnBu",
        mask=enrichment_matrix < threshold,
        ax=ax,
        cbar_kws={'label': '-log10(p-value)'},
        linewidths=0.5,
        linecolor='gray'
    )
    ax.set_title("Cell Type Enrichment (-log10 p-value)")
    ax.set_xlabel("Module")
    ax.set_ylabel("Cell Type")

    return enrichment_matrix, fig





def run_infomap(adata, markov_time=1.0, key_added="infomap"):
    """
    Run Infomap clustering on a Scanpy AnnData object and store results in adata.obs.

    Parameters:
        adata (AnnData): Preprocessed Scanpy object with neighbors computed
        markov_time (float): Infomap Markov time parameter
        key_added (str): Column name to store the clustering labels in adata.obs

    Returns:
        pd.Series: Cluster labels stored in adata.obs[key_added]
    """
    # Convert to NetworkX graph
    try:
        G = nx.from_scipy_sparse_array(adata.obsp['connectivities'])
    except AttributeError:
        G = nx.from_scipy_sparse_matrix(adata.obsp['connectivities'])

    # Build edge list for Infomap
    edges = list(G.edges(data=True))
    src = [e[0] for e in edges]
    trg = [e[1] for e in edges]
    weight = [e[2].get('weight', 1.0) for e in edges]

    n_nodes = adata.n_obs
    im = Infomap(f"--two-level --directed --markov-time {markov_time}")

    for i in range(len(src)):
        im.add_link(int(src[i]), int(trg[i]), weight[i])

    im.run()

    # Extract clustering
    cids = np.zeros(n_nodes)
    for node in im.tree:
        if node.is_leaf:
            cids[node.node_id] = node.module_id

    # Normalize to sequential IDs
    labels = pd.Series(np.unique(cids, return_inverse=True)[1], index=adata.obs_names)

    # Store in adata.obs
    adata.obs[key_added] = labels.astype(str)

    return adata.obs[key_added]