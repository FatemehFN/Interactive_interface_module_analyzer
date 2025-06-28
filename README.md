 # <font color="purple">  READ ME </font>

There are eight functions in `backend.py`. This document will explain the purpose of each function. 

## 1. <font color="lightblue">build_knn_graph(expr_matrix, k)</font>

This function constructs a **k-nearest neighbors graph** from single-cell gene expression data using the [Scanpy](https://scanpy.readthedocs.io/) toolkit. 

#### Parameters:
- `expr_matrix` (`array-like`): A 2D matrix of shape `(n_cells, n_genes)` containing gene expression values.
- `k` (`int`): The number of nearest neighbors to use when building the graph.

#### **Returns:**
- `adata` (`AnnData`): A Scanpy `AnnData` object containing:
  - The original expression matrix
  - PCA-transformed data (`adata.obsm['X_pca']`)
  - kNN graph (`adata.obsp['connectivities']`)
  - UMAP embedding (`adata.obsm['X_umap']`)
  - t-SNE embedding (`adata.obsm['X_tsne']`)

## 2.  <font color="lightblue">run_leiden(adata, resolution)</font>


This function performs **Leiden clustering** on a single-cell dataset using the [Scanpy](https://scanpy.readthedocs.io/) toolkit. 

#### Parameters:

- **adata** (`AnnData`): Annotated data matrix containing neighborhood graph.

- **resolution** (`float`): Determines the granularity of clustering. Higher values yield more clusters.

#### Returns:

- **pandas.Series**: Cluster labels assigned to each cell (`adata.obs['leiden']`), with each label representing a distinct community.

## 3.  <font color="lightblue">compute_eigengenes(expr_matrix, cluster_labels)</font>


This function computes module eigengenes using PCA. 

#### Parameters:

- **expr_matrix** (`pd.DataFrame`):  
  Gene expression matrix with genes as rows and samples as columns.

- **cluster_labels** (`list` or `pd.Series`):  
  Cluster labels assigning each gene to a module. 

#### Returns:

- **`pd.DataFrame`**:  
  Eigengene matrix with samples as rows and modules as columns (named as `Module_0`, `Module_1`, etc.).


## 4.  <font color="lightblue">correlate_eigengenes(eigengenes, phenotypes)</font>


This function computes Pearson correlation coefficients and associated p-values between module eigengenes and phenotype traits.


#### Parameters:

- **eigengenes** (`pd.DataFrame`):  
  A matrix of shape *(samples × modules)*, where each column is a module eigengene.

- **phenotypes** (`pd.DataFrame`):  
  A matrix of shape *(samples × traits)*, where each column represents a phenotype.


#### Returns:

- **`(pd.DataFrame, pd.DataFrame)`**:  
  1. **Correlation matrix**  
  2. **P-value matrix** 



## 5.  <font color="lightblue">correlate_eigengenes_Biweight_midcorrelation(eigengenes, phenotypes)</font>


This function computes biweight midcorrelation and corresponding p-values between module eigengenes and phenotype traits. This function and the previous one both calculate correlations using different methods; running both increases confidence in the robustness of the results.


#### Parameters:

- **eigengenes** (`pd.DataFrame`):  
  A matrix of shape (samples × modules), where each column is a module eigengene.

- **phenotypes** (`pd.DataFrame`):  
  A matrix of shape (samples × traits), where each column is a phenotype.



#### Returns:

- **`(pd.DataFrame, pd.DataFrame)`**:  
  - Correlation matrix  
  - P-value matrix  
 

## 6.  <font color="lightblue">analyze_best_k(expr_matrix, k_values, resolution=1.0)</font>

This function evaluates how varying the number of nearest neighbors (`k`) affects Leiden clustering results.

#### Parameters:

- **`adata`** (`AnnData`):  
 Annotated data matrix containing neighborhood graph.

- **`k_values`** (`list` of `int`):  
  List of `k` values (number of neighbors) to test.

- **`resolution`** (`float`):  
  Resolution parameter for Leiden clustering controlling cluster granularity.


#### Returns:

- Generates and saves plots showing the summary statistics for each tested `k` value.

## 7.  <font color="lightblue">perform_cell_type_enrichment_and_heatmap(expr_df, clusters, threshold=1e-3)</font>

This function utilizes a curated set of brain cell type marker genes (e.g., neurons) to perform enrichment analysis on gene modules, then visualizes the significance in a heatmap.


#### Parameters:

- **`expr_df`** (`pd.DataFrame`):  
  Gene expression DataFrame indexed by gene names (rows) with expression values.

- **`clusters`** (`array-like` or `pd.Series`):  
  Cluster/module assignments for each gene, indexed by gene names.

- **`threshold`** (`float`, optional, default=1e-3):  
  Minimum enrichment score (-log10 p-value) threshold for displaying values on the heatmap. Scores below this are masked.


#### Returns:

- **`enrichment_matrix`** (`pd.DataFrame`):  
  Matrix of enrichment scores (-log10 p-values) for each cell type (rows) across gene modules (columns).

- **`fig`** (`matplotlib.figure.Figure`):  
  Heatmap figure visualizing the enrichment scores with annotations indicating statistical significance.

## 8.  <font color="lightblue">run_infomap(adata, markov_time=1.0, key_added="infomap")</font>


This function performs Infomap community detection on the kNN graph and stores the resulting cluster labels in the AnnData object. It converts the kNN graph to a NetworkX graph and builds an edge list with weights.


#### Parameters:

- **`adata`** (`AnnData`):  

- *`markov_time`** (`float`, default = `1.0`):  
  The Markov time parameter for Infomap, controlling the granularity of the clustering. Higher values may produce more clusters.

- **`key_added`** (`str`, default = `"infomap"`):  
  The name of the column in `adata.obs` where the resulting cluster labels will be stored.



#### Returns:

- **`pd.Series`**:  
  Cluster labels indexed by cell, stored in `adata.obs[key_added]`.


