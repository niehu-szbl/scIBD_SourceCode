import os
import pandas as pd
import numpy as np
import scanpy as sc

## read data
adata = sc.read_h5ad("all_CD8.h5ad")

## special gene clusters
IG = pd.read_table("/home/niehu/niehu/database/00.genes/excluded/ig.txt", header = None)
TCR = pd.read_table("/home/niehu/niehu/database/00.genes/excluded/tr3.txt", header = None)
CC = pd.read_table("/home/niehu/niehu/database/00.genes/excluded/regev_lab_cell_cycle_genes.txt", header=None)
IG = IG.loc[:,0].to_list()
TCR = TCR.loc[:,0].to_list()
CC = CC.loc[:,0].to_list()
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rp'] = adata.var_names.str.match('^RP[LS]') | adata.var_names.str.match('^RP[0-9]+-')
mt = adata.var.index[ adata.var.mt].to_list()
rp = adata.var.index[ adata.var.rp].to_list()
exclude = IG + TCR + mt + rp

## filter cells
#pd.DataFrame( adata.obs.sampleName.value_counts()).to_csv("cell_per_sample.csv", sep = "\t", header=None)
select = adata.obs.sampleName.value_counts() > 10
select = select.index[select].to_list()
adata = adata[ adata.obs.sampleName.isin(select)]

## filter genes
select_gene = [gene not in exclude for gene in adata.var.index]
adata = adata[ :,select_gene]

## run scanpy
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## find hvg
select = adata.obs.sampleName.value_counts() > 10
select = select.index[select].to_list()
adata_tmp = adata[ adata.obs.sampleName.isin(select)]
sc.pp.highly_variable_genes(adata_tmp, n_top_genes = 2000, batch_key = "sampleName")
#pd.Series(adata_tmp.var_names[adata_tmp.var.highly_variable].to_list()).to_csv("hvg_n100.csv", header=None, index=False)
hvg = adata_tmp.var_names[adata_tmp.var.highly_variable].to_list()
del adata_tmp

## find and assign hvg
sc.pp.highly_variable_genes(adata, n_top_genes = 2000, batch_key = "sampleName")
adata.var.highly_variable = adata.var.index.isin(hvg)
#pd.Series(adata.var_names[adata.var.highly_variable].to_list()).to_csv("hvg_n10.csv", header=None, index=False)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
adata.write_h5ad("all_CD8.n10.h5ad")
