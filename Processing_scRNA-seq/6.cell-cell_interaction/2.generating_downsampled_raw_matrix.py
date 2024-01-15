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
#import scvelo as scv
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/total_merge_filtered.h5ad')
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')
meta_data = pd.read_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/subset_total_500_meta.csv')
meta_data = meta_data[['barcode','celltype']]
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata = adata[adata.obs['barcode'].isin(meta_data['barcode']), :]
adata.obs = pd.merge(adata.obs, meta_data, on = 'barcode')
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
