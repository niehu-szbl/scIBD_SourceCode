import bbknn
import scanpy as sc

adata = sc.read_h5ad("Endothelial.bySample.n10.h5ad")

bbknn.bbknn(adata, batch_key='sampleName', neighbors_within_batch=3, metric='euclidean', n_pcs=30)
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution = 0.3)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write('Endothelial.bySample.n10.bbknn_umap_leiden_marker.res0.3.h5ad')

sc.tl.leiden(adata, resolution = 0.6)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write('Endothelial.bySample.n10.bbknn_umap_leiden_marker.res0.6.h5ad')

sc.tl.leiden(adata, resolution = 0.9)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write('Endothelial.bySample.n10.bbknn_umap_leiden_marker.res0.9.h5ad')

sc.tl.leiden(adata, resolution = 1.2)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write('Endothelial.bySample.n10.bbknn_umap_leiden_marker.res1.2.h5ad')

sc.tl.leiden(adata, resolution = 1.5)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
adata.write('Endothelial.bySample.n10.bbknn_umap_leiden_marker.res1.5.h5ad')
