import bbknn
import scanpy as sc

adata = sc.read_h5ad("ibd_20220111.h5ad")
bbknn.bbknn(adata, batch_key='sampleName', neighbors_within_batch=3, metric='euclidean', n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.1)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write('ibd_20220111.bbknn_umap_leiden_marker.h5ad')
