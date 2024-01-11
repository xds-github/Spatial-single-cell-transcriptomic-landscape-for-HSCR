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
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
scv.logging.print_version()
ID="J20060200-5"
ID2= "CT1"
group = "CT"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20070130-5"
ID2= "CT2"
group = "CT"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20100206-5"
ID2= "CT7"
group = "CT"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21010081-5"
ID2= "CT10"
group = "CT"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21030056-5"
ID2= "CT11"
group = "CT"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21030122-5"
ID2= "CT12"
group = "CT"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20080085-5"
ID2= "HL25K"
group = "HLG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20080086-5"
ID2= "HL25X"
group = "HLA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21030052-5"
ID2= "HL41K"
group = "HLG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21030051-5"
ID2= "HL41X"
group = "HLA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21040036-5"
ID2= "HL45K"
group = "HLG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21040035-5"
ID2= "HL45X"
group = "HLA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21050062-5"
ID2= "HL55K"
group = "HLG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J21050061-5"
ID2= "HL55X"
group = "HLA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20070149-5"
ID2= "HS14K"
group = "HSG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20070150-5"
ID2= "HS14X"
group = "HSA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20090072-5"
ID2= "HS26K"
group = "HSG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20090073-5"
ID2= "HS26X"
group = "HSA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20090107-5"
ID2= "HS28K"
group = "HSG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20090106-5"
ID2= "HS28X"
group = "HSA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20100208-5"
ID2= "HS33K"
group = "HSG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20100207-5"
ID2= "HS33X"
group = "HSA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20110101-5"
ID2= "HS35K"
group = "HSG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20110100-5"
ID2= "HS35X"
group = "HSA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20110165-5"
ID2= "HS36K"
group = "HSG"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

ID="J20110164-5"
ID2= "HS36X"
group = "HSA"
adata1 = scv.read('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/velocyto/'+ID+'.loom', cache=True)
results_file = '/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/'+ID2+'.h5ad'  # the file that will store the analysis results
adata2 = sc.read_10x_mtx('/share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/'+ID+'/outs/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=True)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.obs['barcode'] = adata2.obs._stat_axis.values.tolist()
adata = scv.utils.merge(adata2, adata1)
adata.obs.rename(index=adata.obs['barcode'],inplace=True)
adata.obs = adata.obs.drop('barcode', axis=1)
print(adata)
adata.obs["ID"] = ID2
adata.obs["group"] = group
adata.write(results_file)

CT1 = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/CT1.h5ad")
CT2 = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/CT2.h5ad")
CT7 = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/CT7.h5ad")
CT10 = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/CT10.h5ad")
CT11 = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/CT11.h5ad")
CT12 = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/CT12.h5ad")
HL25K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL25K.h5ad")
HL25X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL25X.h5ad")
HL41K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL41K.h5ad")
HL41X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL41X.h5ad")
HL45K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL45K.h5ad")
HL45X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL45X.h5ad")
HL55K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL55K.h5ad")
HL55X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HL55X.h5ad")
HS14K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS14K.h5ad")
HS14X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS14X.h5ad")
HS26K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS26K.h5ad")
HS26X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS26X.h5ad")
HS28K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS28K.h5ad")
HS28X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS28X.h5ad")
HS33K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS33K.h5ad")
HS33X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS33X.h5ad")
HS35K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS35K.h5ad")
HS35X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS35X.h5ad")
HS36K = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS36K.h5ad")
HS36X = sc.read("/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/HS36X.h5ad")

adata = CT1.concatenate(CT2,CT7,CT10,CT11,CT12,HL25K,HL25X,HL41K,HL41X,HL45K,HL45X,HL55K,HL55X,HS14K,HS14X,HS26K,HS26X,HS28K,HS28X,HS33K,HS33X,HS35K,HS35X,HS36K,HS36X)
del CT1,CT2,CT7,CT10,CT11,CT12,HL25K,HL25X,HL41K,HL41X,HL45K,HL45X,HL55K,HL55X,HS14K,HS14X,HS26K,HS26X,HS28K,HS28X,HS33K,HS33X,HS35K,HS35X,HS36K,HS36X
#adata1 = adata.raw.to_adata()
#adata = adata.raw.to_adata()
sc.pp.filter_genes(adata, min_cells=3)
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')

adata=sc.read('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/raw_merge.h5ad')
#adata = adata.raw.to_adata()
#sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 25, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=[ 'group'])
sc.external.pp.bbknn(adata, batch_key='ID')
sc.tl.umap(adata)
sc.pl.umap(adata, color=[ 'group'])
sc.tl.leiden(adata,key_added="clusters")
rcParams['figure.figsize'] = 5,5
sc.pl.umap(adata, color=[ 'clusters'], legend_loc= 'on data')
adata.write('/share/home/xudeshu/scanpy_dic/HSCR/h5ad_out/total_merge.h5ad')

sc.pl.umap(adata, color=[ 'clusters'], legend_loc= 'on data')

sc.pl.dotplot(adata, ["PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1","CLEC9A","CD1C",'LAMP3',"CD69","CD38","CXCR3","EPCAM","PDPN","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",'S100A8','S100A9','PTPRZ1'], groupby='clusters',cmap ="Spectral_r")
sc.pl.dotplot(adata, ['PECAM1','LYVE1',"PTPRC","CD3E",'CD3D',"CD8A","CD4","KIT","CD19","MME","IGHA1","SDC1","CD14","S100B","GPM6B","TPSAB1",
                      "CLEC9A","CD1C",'LAMP3',"CD69","CD38","CXCR3","EPCAM","PDPN","COL3A1","CXCL14","CFD","LIF","THBS4","PECAM1","BCAM","MYH11",
                      'S100A8','S100A9','PTPRZ1','NOTCH3','RGS5','CLEC4C'], groupby='clusters',cmap ="Spectral_r")
# -------------------------------------------------------- subset by cell lineages --------------------------------------------------------------
# Epi
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['14','21','30']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Epi_meta_data.csv", sep=',')

# FLC
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['1','20','22','3','9','19']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/FLC_meta_data.csv", sep=',')

# Glial
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['24']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Glial_meta_data.csv", sep=',')

# Myloid
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['10','13','25','27','29']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Myloid_meta_data.csv", sep=',')

# T and ILC cell
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['6','7','16','18','23','29','30']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/T_ILC_meta_data.csv", sep=',')

# B cell
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['0','2','4','5','12','17','30','8','11']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/B_meta_data.csv", sep=',')

# subset Pericyte
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['15']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Peri_meta_data.csv", sep=',')

# subset Endo
adata.obs['barcode'] = adata.obs._stat_axis.values.tolist()
adata1 = adata[adata.obs['clusters'].isin(['26','28']), :]
meta_data =adata1.obs
meta_data[['barcode']].to_csv( "/share/home/xudeshu/scanpy_dic/HSCR/temp_meta/Endo_meta_data.csv", sep=',')
