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
#adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/total_merge_filtered.h5ad')
# -------------------------------------------------------------- filtering total -----------------------------------------------------------------
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/total_merge.h5ad')
meta_data = pd.read_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/total_meta_data.csv') # Read in a merged filtered meta data
meta_data = meta_data[['barcode','celltype']]
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata = adata[adata.obs['barcode'].isin(meta_data['barcode']), :]
adata.obs = pd.merge(adata.obs, meta_data, on = 'barcode')
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
sc.pl.umap(adata, color=[ 'celltype'])
# Set celltypes order
celltype_order = ['Telocyte_VSTM2A','Telocyte_NPY','FLC_ACTA2','FLC_SCARA5','FLC_CCL19','FLC_GREM1','FLC_APOD','FLC_CLDN1','FLC_KCNN3',
                  'Pericyte_RGS5','Pericyte_MYH11','Pericyte_CCL2','EC_Lym','EC_Venous','EC_Arterial',
                  'Epi_stem', 'TA','Enterocytes', 'Enterocytes_SLC26A3','Enterocytes_LYZ', 'Enterocytes_BEST4', 'Goblet_RETNLB', 'Goblet_ACHE',  'Tuft', 'Enteroendocrine',
                  'Glial',
                 'GCB_MKI67','GCB','ProB','Brm_IGL','Brm_FCRL4','Brm_CD267','Brm','Plasmablast','Plasma',
                  'LTi','Th17','Th2','Treg_like','Tfh_like','CD4T_naive','CD8T_naive','CD8Tcm','CD8Tgd_ZNF683','CD8Tgd_CD39','CD8Trm','CD8Tem','CD8Tem_CX3CR1','T_cycling','gdT','NK',
                  'M_IL1B','M_C1Q','M_LYVE1','M_APOE','M_cycling','cDC2','cDC1','DC_LAMP3','pDC','Mast']
adata.obs['celltype'] = adata.obs['celltype'].cat.set_categories(celltype_order)
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/total_merge_filtered.h5ad')

adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
meta_data =adata.obs
meta_data[['barcode','celltype','ID','group']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/total_meta_data2.csv", sep=',')
# ------------------------------------------------------- DEG analysis -------------------------------------------------------------------------
# Find markers of the clusters
sc.tl.rank_genes_groups(adata, 'celltype', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
result = adata.uns['rank_genes_groups']
sc.pl.umap(adata, color=[ 'celltype'], legend_loc = "on data", palette = "Set3_r")
groups = result['names'].dtype.names
total_markers = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj','logfoldchanges','scores']})
total_markers.to_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/sc_all_markers.csv')
total_markers.head(10)

adata_temp = adata[adata.obs['celltype'].isin(['Pericyte_RGS5','Pericyte_MYH11','Pericyte_CCL2']), :]
#adata_temp = adata
sc.tl.rank_genes_groups(adata_temp, 'celltype', groups=['Pericyte_CCL2'], method='wilcoxon')
sc.pl.rank_genes_groups(adata_temp, n_genes=25, sharey=False)
result = adata_temp.uns['rank_genes_groups']
groups = result['names'].dtype.names
total_markers = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj','logfoldchanges','scores']})
total_markers.to_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/Pericyte_CCL2_v_restPer.csv')
total_markers.head(10)

adata_temp = adata[adata.obs['celltype'].isin(['FLC_APOD']), :]
sc.tl.rank_genes_groups(adata_temp, 'group2', method='wilcoxon')
sc.pl.rank_genes_groups(adata_temp, n_genes=25, sharey=False)
result = adata_temp.uns['rank_genes_groups']
sc.pl.umap(adata_temp, color=[ 'celltype'], legend_loc = "on data", palette = "Set3_r")
groups = result['names'].dtype.names
total_markers = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj','logfoldchanges','scores']})
total_markers.to_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/sc_all_markers.csv')
total_markers.head(10)


