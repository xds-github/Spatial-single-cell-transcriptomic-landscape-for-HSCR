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
#adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/FCL_out/FLC_processed.h5ad')
#scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
#scv.settings.presenter_view = True  # set max width size for presenter view
#scv.settings.set_figure_params('scvelo')  # for beautified visualization
pair11 = ["#C294B9","#E2D5E6","#A60B75","#6A371B","#B17D30","#2E796D","#B4B5B5","#C37284","#1E2963","#D6A128","#ECCB20"]
#adata = sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/T_ILC_processed-v2.h5ad')
df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/T_ILC_meta_data.csv")
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
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/T_ILC.h5ad')

adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/T_ILC.h5ad')

sc.pl.dotplot(adata, ['PECAM1','LYVE1',"PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1",
                      "CLEC9A","CD1C",'LAMP3',"CD69","CD38","CXCR3","EPCAM","PDPN","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",
                      'S100A8','S100A9','NOTCH3','RGS5','CD14','MKI67'], groupby='clusters',cmap ="Spectral_r")

sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data")

marker_genes = {'Marker genes':['KIT','IL23R','RORC','CCR6','PTGDR2','IL13','IL17RB','IL2RA','FOXP3','CTLA4','PDCD1','CD4','CD3D','CD3E','KLRB1','SELL','CCR7','CD8A','CD8B','ZNF683','TRDV1','ENTPD1','KIR2DL4','ITGAE',
                'ITGA1','CD160','GNLY','PRF1',
                'CX3CR1','FCGR3A','MKI67','TRDV2','KLRF1','CMC1','KLRC1','EOMES','XCL1','XCL2'],
                'Hormone receptors':['ADRB2','P2RX1','P2RX4','P2RX7','CHRNB1','SPPL3'],
               'Granzyme':['GZMK','GZMA','GZMB','GZMH','GZMM'],
                'MHC Class II':['HLA-DRB1','HLA-DPB1','HLA-DPA1','HLA-DRB5','HLA-DQB1','HLA-DQA1','HLA-DMA','SLC4A10']}
#sc.pl.dotplot(adata, marker_genes, groupby='celltype', cmap='YlGnBu_r',save='T_marker.pdf')
sc.pl.dotplot(adata, marker_genes, groupby='clusters', cmap='YlGnBu_r')

# subset Endo
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['0','1','2','3','4','5','6','8','9','10','11','12','15','17']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/T_ILC_meta_data_v2.csv", sep=',')

df_news = pd.read_csv("/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/T_ILC_meta_data_v2.csv")
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
sc.tl.leiden(adata,key_added="clusters", resolution = 1.5)
sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data")
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/T_ILC_v2.h5ad')

adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/temp/T_ILC_v2.h5ad')

sc.pl.dotplot(adata, ['PECAM1','LYVE1',"PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1",
                      "CLEC9A","CD1C",'LAMP3',"CD69","CD38","CXCR3","EPCAM","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",
                      'S100A8','S100A9','NOTCH3','RGS5','CD14','MKI67'], groupby='clusters',cmap ="Spectral_r")

#sc.tl.leiden(adata,key_added="clusters", resolution = 1.5)
sc.pl.umap(adata, color=[ 'clusters'], legend_loc = "on data")
sc.pl.umap(adata, color=[ 'clusters'])

marker_genes = {'Marker genes':['KIT','IL23R','RORC','CCR6','PTGDR2','IL13','IL17RB','IL2RA','FOXP3','CTLA4','PDCD1','CD4','CD3D','CD3E','KLRB1','SELL','CCR7','CD8A','CD8B','ZNF683','TRDV1','ENTPD1','KIR2DL4','ITGAE',
                'ITGA1','CD160','GNLY','PRF1',
                'CX3CR1','FCGR3A','MKI67','TRDV2','KLRF1','CMC1','KLRC1','EOMES','XCL1','XCL2'],
                'Hormone receptors':['ADRB2','P2RX1','P2RX4','P2RX7','CHRNB1','SPPL3'],
               'Granzyme':['GZMK','GZMA','GZMB','GZMH','GZMM'],
                'MHC Class II':['HLA-DRB1','HLA-DPB1','HLA-DPA1','HLA-DRB5','HLA-DQB1','HLA-DQA1','HLA-DMA','SLC4A10']}
#sc.pl.dotplot(adata, marker_genes, groupby='celltype', cmap='YlGnBu_r',save='T_marker.pdf')
sc.pl.dotplot(adata, marker_genes, groupby='clusters', cmap='YlGnBu_r')
# Anotate the celltype 
adata.obs['celltype'] = adata.obs['clusters']
categ = adata.obs['celltype'].cat.categories
new_cluster_names = ['CD8Tgd_ZNF683','CD4T_naive','CD8Tem','CD8T_naive','NK','Th17','CD8Tgd_ZNF683','CD8Tgd_CD39','LTi','Treg_like','CD8Trm','Tfh_like','CD8Tem_CX3CR1','T_cycling','gdT','CD8Trm','Th2','CD8Tcm','CD8Tcm','Polluted']
adata.obs['celltype'] = adata.obs['celltype'].cat.add_categories(np.unique(new_cluster_names))
for i in range(len(adata.obs['clusters'].cat.categories)):
    adata.obs['celltype'][adata.obs['clusters']==adata.obs['clusters'].cat.categories[i]] = new_cluster_names[i]
adata.obs['celltype'] = adata.obs['celltype'].cat.remove_categories(categ)
del categ

# Set celltypes order
celltype_order = ['LTi','Th17','Th2','Treg_like','Tfh_like','CD4T_naive','CD8T_naive','CD8Tcm','CD8Tgd_ZNF683','CD8Tgd_CD39','CD8Trm','CD8Tem','CD8Tem_CX3CR1','T_cycling','gdT','NK','Polluted']
adata.obs['celltype'] = adata.obs['celltype'].cat.set_categories(celltype_order)
adata = adata[adata.obs['celltype'].isin(['LTi','Th17','Th2','Treg_like','Tfh_like','CD4T_naive','CD8T_naive','CD8Tcm','CD8Tgd_ZNF683','CD8Tgd_CD39','CD8Trm','CD8Tem','CD8Tem_CX3CR1','T_cycling','gdT','NK']), :]
marker_genes = {'Marker genes':['KIT','IL23R','RORC','CCR6','PTGDR2','IL13','IL17RB','IL2RA','FOXP3','CTLA4','PDCD1','CD4','CD3D','CD3E','KLRB1','SELL','CCR7','CD8A','CD8B','ZNF683','TRDV1','ENTPD1','KIR2DL4','ITGAE',
                'ITGA1','CD160','GNLY','PRF1',
                'CX3CR1','FCGR3A','MKI67','TRDV2','KLRF1','CMC1','KLRC1','EOMES','XCL1','XCL2'],
                'Hormone receptors':['ADRB2','P2RX1','P2RX4','P2RX7','CHRNB1','SPPL3'],
               'Granzyme':['GZMK','GZMA','GZMB','GZMH','GZMM'],
                'MHC Class II':['HLA-DRB1','HLA-DPB1','HLA-DPA1','HLA-DRB5','HLA-DQB1','HLA-DQA1','HLA-DMA','SLC4A10']}
#sc.pl.dotplot(adata, marker_genes, groupby='celltype', cmap='YlGnBu_r',save='T_marker.pdf')
sc.pl.dotplot(adata, marker_genes, groupby='celltype', cmap='YlGnBu_r')
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/h5ad_file/T_ILC_processed.h5ad')
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
meta_data =adata.obs
meta_data[['barcode','celltype']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/T_ILC_meta_data.csv", sep=',')
