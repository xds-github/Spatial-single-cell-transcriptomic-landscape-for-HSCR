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
#adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/Endo_peri_processed.h5ad')
#scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
#scv.settings.presenter_view = True  # set max width size for presenter view
#scv.settings.set_figure_params('scvelo')  # for beautified visualization
pair11 = ["#FFED6F","#BC80BD","#FCCDE5","#8DD3C7","#BEBADA","#80B1D3","#B4B5B5","#C37284","#1E2963","#D6A128","#ECCB20"]
df_temp1 = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Peri_meta_data.csv")
df_temp2 = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Endo_meta_data.csv")
df_news = pd.concat([df_temp1,df_temp2],axis=0,ignore_index=True)
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
sc.pl.umap(adata, color=['group'])
# Exclude the batch effect
sc.external.pp.bbknn(adata, batch_key='ID')
# Rerun the UMAP: no index need to be changed
sc.tl.umap(adata)
sc.pl.umap(adata, color=[ 'group'])
# Find clusters
sc.tl.leiden(adata,key_added="clusters", resolution = 1)
sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data")
sc.pl.dotplot(adata, ['PECAM1','LYVE1',"PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1",
                      "CLEC9A","CD1C","CD69","CD38","CXCR3","EPCAM","PDPN","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",
                      'S100A8','S100A9','NOTCH3','RGS5','PDGFRA','PDGFRB','RGS5','NOTCH3','PDGFRB','CSPG4'], groupby='clusters',cmap ="Spectral_r")
# subset Endo
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['0','1','2','3','4','5','6','7','14']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Peri_Endo_meta_data.csv", sep=',')
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Peri_Endo_meta_data.csv")
adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/Endo_peri')
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
sc.pl.umap(adata, color=['group'])
# Exclude the batch effect
sc.external.pp.bbknn(adata, batch_key='ID')
# Rerun the UMAP: no index need to be changed
sc.tl.umap(adata)
sc.pl.umap(adata, color=[ 'group'])
# Find clusters
sc.tl.leiden(adata,key_added="clusters", resolution = 1)
sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data")
# Anotate the celltype 
adata.obs['celltype'] = adata.obs['clusters']
categ = adata.obs['celltype'].cat.categories
new_cluster_names = ['Pericyte_RGS5','Pericyte_MYH11','Pericyte_CCL2','Pericyte_RGS5','EC_Lym','EC_Arterial','EC_Venous','EC_Arterial']
adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(np.unique(new_cluster_names))
for i in range(len(adata.obs['clusters'].cat.categories)):
    adata.obs['celltype'][adata.obs['clusters']==adata.obs['clusters'].cat.categories[i]] = new_cluster_names[i]
adata.obs['celltype'] = adata.obs['celltype'].cat.remove_categories(categ)
del categ
# Set celltypes order
celltype_order = ['Pericyte_MYH11','Pericyte_CCL2','Pericyte_RGS5','EC_Lym','EC_Venous','EC_Arterial']
adata.obs['celltype'] = adata.obs['celltype'].cat.set_categories(celltype_order)
sc.pl.umap(adata, color=[ 'celltype'], legend_loc = "on data", palette =['#B2D2CE','#D58889','#7594A7','#FCCDE5','#94D5CA','#C18AC2'], save='vessel_UMAP.pdf')
marker_genes = ['NOTCH3','ACTA2','RGS5','PDGFRB','MYH11','ACTG2','MUSTN1','MRGPRF','CCL2','PROCR','PECAM1','LYVE1','PROX1','CCL14','ACKR1','APLNR','RGCC','HEY1','CXCL12','TM4SF1','IGFBP5','ANGPTL1','ANGPTL2',
               'VEGFA','VEGFB','VEGFC','PDGFD','NOTCH1','NOTCH2','NOTCH3','ACTG2','RGS5','GJA4','GJC1']
marker_genes = ['MYH11','ACTG2','MUSTN1','NOTCH3','ACTA2','RGS5','PDGFRB','CCL2','PROCR','PECAM1','LYVE1','PROX1','CCL14','ACKR1','APLNR','RGCC','HEY1','CXCL12']
ax = sc.pl.stacked_violin(adata, marker_genes, groupby='celltype', swap_axes=False, dendrogram=False,categories_order = None, rotation=90,row_palette= ['#B2D2CE','#7594A7','#D58889','#FCCDE5','#94D5CA','#C18AC2'],save='pericyte_endo_markers.pdf')
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/Endo_peri_processed.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
meta_data =adata.obs
meta_data[['barcode','celltype','group','ID']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/Endo_peri_meta_data.csv", sep=',')
adata1 = sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/Endo_peri_meta_data.csv")
adata1.obs['barcode'] = adata1.obs._stat_axis.values.tolist()
adata1 = adata1[adata1.obs['barcode'].isin(df_news['barcode']), :]
adata1.obs = adata1.obs.drop('barcode', axis=1)
#adata = adata.raw.to_adata()
mat=pd.DataFrame(data=adata1.X.todense(),index=adata1.obs_names,columns=adata1.var_names)
mat.to_csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/Endo_peri_raw_matrix_230921.csv')
