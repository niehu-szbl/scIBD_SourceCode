# Prepare gene expression matrix for each selected study

# Combine scRNA-seq datasets, and run PCA on local computer
# run integrate.ipynb
# output ibd_20220111.h5ad

# Integrate these datasets on a computer cluster
# set resolution = 0.1 when run leiden clustring
# run integrate_all.py
# output ibd_20220111.bbknn_umap_leiden_marker.h5ad

# Round 1 clustering
# cluster all cells into 
# run round1_clustering.ipynb
# output ibd_20220111.bbknn_umap_leiden_marker.adata.obs.csv

# Then do the second round clustering
# cd round2_clustering/Neural, integrate neural cells
# cd round2_clustering/Endothelial, integrate endothelial cells
# cd round2_clustering/Epithelial, integrate epithelial cells
# cd round2_clustering/Mesenchymal, integrate mesenchymal cells
# cd round2_clustering/B_plasma, integrate B cells and plasma cells
# cd round2_clustering/Myeloid, integrate myeloid cells
# cd round2_clustering/B_plasma, integrate B cells and plasma cells
# cd round2_clustering/T_NK, integrate T/NK cells