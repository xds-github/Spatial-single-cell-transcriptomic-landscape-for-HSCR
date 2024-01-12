import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import loompy as lp
import matplotlib.pyplot as pl
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
import seaborn as sns
#import csv
from pandas.core.frame import DataFrame
#import scvelo as scv
#adata = sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/FLC_processed-v2.h5ad')
#scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
#scv.settings.presenter_view = True  # set max width size for presenter view
#scv.settings.set_figure_params('scvelo')  # for beautified visualization
pair11 = ["#C294B9","#E2D5E6","#A60B75","#6A371B","#B17D30","#2E796D","#B4B5B5","#C37284","#1E2963","#D6A128","#ECCB20"]
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/FLC_meta_data.csv")
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
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/FLC.h5ad')

adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/FLC.h5ad')

sc.pl.dotplot(adata, ['PECAM1','LYVE1',"PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1",
                      "CLEC9A","CD1C",'LAMP3',"CD69","CD38","CXCR3","EPCAM","PDPN","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",
                      'S100A8','S100A9','NOTCH3','RGS5'], groupby='clusters',cmap ="Spectral_r")

sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data")

# subset Endo
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['0','1','2','3','4','5','6','8','9','10','11','13']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/FLC_meta_data_v2.csv", sep=',')

df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/FLC_meta_data_v2.csv")
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
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/FLC_v2.h5ad')

adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/FLC_v2.h5ad')

sc.pl.dotplot(adata, ['PECAM1','LYVE1',"PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1",
                      "CLEC9A","CD1C",'LAMP3',"CD69","CD38","CXCR3","EPCAM","PDPN","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",
                      'S100A8','S100A9','NOTCH3','RGS5'], groupby='clusters',cmap ="Spectral_r")

marker_genes = ['NOTCH3','RGS5','NDUFA4L2','ACTA2','MYH11','ZFHX3','POSTN','PDGFRA','VSTM2A','NRG1','NPY','CD81','APOE','SOSTDC1','HHIP','ADAMDEC1','SCARA5','CCL19','CD74','RBP5','C3','GPNMB','GREM1','SFRP2','CPZ','MFAP5','BCHE','KCNN3','THBS4','SCN7A','APOD','PI16','CXCL12']
sc.pl.dotplot(adata, marker_genes, groupby='clusters')

# Find markers of the clusters
sc.tl.rank_genes_groups(adata, 'clusters', method='FLC_raw_downsample.h5seurat')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
result = adata.uns['rank_genes_groups']
sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data", palette = "Set3_r")
groups = result['names'].dtype.names
total_markers = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals_adj','logfoldchanges']})
#total_markers.to_csv('/data/scanpy/Epi_markers.csv')
total_markers.head(10)
# Anotate the celltype 
adata.obs['celltype'] = adata.obs['clusters']
categ = adata.obs['celltype'].cat.categories
new_cluster_names = ['FLC_GREM1','FLC_SCARA5','FLC_APOD','FLC_SCARA5','Telocyte_VSTM2A','FLC_GREM1','FLC_KCNN3','FLC_ACTA2','Telocyte_NPY','FLC_GREM1','FLC_CLDN1','FLC_APOD','FLC_CCL19']
adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(np.unique(new_cluster_names))
for i in range(len(adata.obs['clusters'].cat.categories)):
    adata.obs['celltype'][adata.obs['clusters']==adata.obs['clusters'].cat.categories[i]] = new_cluster_names[i]
adata.obs['celltype'] = adata.obs['celltype'].cat.remove_categories(categ)
del categ
# Set celltypes order
celltype_order = ['Telocyte_VSTM2A','Telocyte_NPY','FLC_ACTA2','FLC_SCARA5','FLC_CCL19','FLC_GREM1','FLC_APOD','FLC_CLDN1','FLC_KCNN3']
adata.obs['celltype'] = adata.obs['celltype'].cat.set_categories(celltype_order)
marker_genes = ['POSTN','PDGFRA','VSTM2A','NRG1','NPY','APOE','ACTA2','MYH11','ADAMDEC1','SCARA5','FABP5','CCL11','CCL19','CD74','RBP5','C3','C7','GREM1','IGF2','PCOLCE2','SFRP2','CPZ','APOD','PI16','CLDN1','CYP1B1','KCNN3','THBS4','SCN7A']
ax = sc.pl.stacked_violin(adata, marker_genes, groupby='celltype', swap_axes=False, dendrogram=False,categories_order = None, rotation=90,
                          row_palette= ["#C7B5D1",'#DE7F20','#ADCBDA','#EABC80','#B3D395','#EFA6A5','#3777A2','#41923B','#CA3335'], save='_FLC_markers.pdf')
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/FLC_processed.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
meta_data =adata.obs
meta_data[['barcode','celltype']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/FLC_meta_data.csv", sep=',')
