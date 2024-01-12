import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import loompy as lp
import matplotlib.pyplot as pl
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
import seaborn as sns
import csv
from pandas.core.frame import DataFrame
#import scvelo as scv
#adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/FCL_out/FLC_processed.h5ad')
#scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
#scv.settings.presenter_view = True  # set max width size for presenter view
#scv.settings.set_figure_params('scvelo')  # for beautified visualization
pair11 = ["#C294B9","#E2D5E6","#A60B75","#6A371B","#B17D30","#2E796D","#B4B5B5","#C37284","#1E2963","#D6A128","#ECCB20"]
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Glial_meta_data.csv")
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata = adata[adata.obs['barcode'].isin(df_news['barcode']), :]
adata.obs = adata.obs.drop('barcode', axis=1)
sc.pp.filter_genes(adata, min_cells=10)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 25, :]
# Normalize and log data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Pick highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
# Regress and scale data
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
# Set the optimal
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
# Find clusters
sc.tl.leiden(adata,key_added="clusters", resolution = 0.5)
# subset Endo
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['0','1','2','3','4','7']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Glial_meta_data.csv_v2", sep=',')
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/Glial.h5ad')
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Glial_meta_data.csv_v2")
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata = adata[adata.obs['barcode'].isin(df_news['barcode']), :]
adata.obs = adata.obs.drop('barcode', axis=1)
sc.pp.filter_genes(adata, min_cells=10)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 25, :]
# Normalize and log data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Pick highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
# Regress and scale data
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
# Set the optimal
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
# subset Endo
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata = adata[adata.obs['clusters'].isin(['0','1','2','3','5','7']), :]
#meta_data =adata1.obs
#meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Glial_meta_data.csv_v3", sep=',')
#adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/Glial.h5ad')
# Find clusters
sc.tl.leiden(adata,key_added="clusters", resolution = 0.5)
sc.pl.umap(adata, color=['group','clusters'])
# Anotate the celltype 
adata.obs['celltype'] = adata.obs['clusters']
categ = adata.obs['celltype'].cat.categories
new_cluster_names = ['Glial','Glial','Glial','Glial','Glial','Glial']
adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(np.unique(new_cluster_names))
for i in range(len(adata.obs['clusters'].cat.categories)):
    adata.obs['celltype'][adata.obs['clusters']==adata.obs['clusters'].cat.categories[i]] = new_cluster_names[i]
adata.obs['celltype'] = adata.obs['celltype'].cat.remove_categories(categ)
del categ
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/glial_processed.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
meta_data =adata.obs
meta_data[['barcode','celltype']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/Glial_meta_data.csv", sep=',')
