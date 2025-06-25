### 🔗 `build_knn_graph(expr_matrix, k)`

This function constructs a **k-nearest neighbors graph** from single-cell gene expression data using the [Scanpy](https://scanpy.readthedocs.io/) toolkit. It also computes low-dimensional embeddings for visualization.

#### **Parameters:**
- `expr_matrix` (`array-like`): A 2D matrix of shape `(n_cells, n_genes)` containing gene expression values.
- `k` (`int`): The number of nearest neighbors to use when building the graph.

#### **Returns:**
- `adata` (`AnnData`): A Scanpy `AnnData` object containing:
  - The original expression matrix
  - PCA-transformed data (`adata.obsm['X_pca']`)
  - kNN graph (`adata.obsp['connectivities']`)
  - UMAP embedding (`adata.obsm['X_umap']`)
  - t-SNE embedding (`adata.obsm['X_tsne']`)

#### **Example Usage:**
```python
adata = build_knn_graph(expr_matrix=my_data, k=15)
sc.pl.umap(adata)
