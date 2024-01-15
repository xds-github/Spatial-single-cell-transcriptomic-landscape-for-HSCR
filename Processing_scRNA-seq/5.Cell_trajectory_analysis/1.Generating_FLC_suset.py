import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
#import loompy as lp
import matplotlib.pyplot as pl
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
import seaborn as sns
import csv
from pandas.core.frame import DataFrame
# -------------------------------------------------------------- extract the raw matirx -----------------------------------------------------------------
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')
meta_data = pd.read_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/total_meta_data2.csv')
meta_data = meta_data[['barcode','celltype']]
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata = adata[adata.obs['barcode'].isin(meta_data['barcode']), :]
adata.obs = pd.merge(adata.obs, meta_data, on = 'barcode')
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
adata = adata[adata.obs['celltype'].isin(['Telocyte_VSTM2A','Telocyte_NPY','FLC_ACTA2','FLC_SCARA5','FLC_CCL19','FLC_GREM1','FLC_APOD','FLC_CLDN1','FLC_KCNN3']), :]
adata = sc.pp.subsample(adata, n_obs = 20000, copy = True)
#sc.pl.umap(adata, color=[ 'celltype'])
# Set celltypes order
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/FLC_raw_downsample.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
meta_data =adata.obs
meta_data[['barcode','celltype']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/FLC_meta_downsample_data.csv", sep=',')
# ------------------------------------------------------------- extract the UMAP ----------------------------------------------------------------------
adata1=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/FLC_processed.h5ad')
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/FLC_meta_downsample_data.csv")
adata1.obs['barcode'] = adata1.obs._stat_axis.values.tolist()
adata1 = adata1[adata1.obs['barcode'].isin(df_news['barcode']), :]
adata1.obs = adata1.obs.drop('barcode', axis=1)
umapDF = pd.DataFrame(adata1.obsm['X_umap'], columns=['_X', '_Y'])
umapDF['barcode'] = adata1.obs._stat_axis.values.tolist()
umapDF.rename(index=umapDF['barcode'],inplace=True)
umapDF = umapDF.drop('barcode', axis=1)
umapDF.to_csv('FLC_UMAP.csv')
