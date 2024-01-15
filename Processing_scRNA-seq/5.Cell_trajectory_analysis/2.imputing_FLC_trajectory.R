library(Seurat)
library(reshape2)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(SeuratDisk)
Convert("J:/scanpy_out/final_ann/monocle3/FLC_raw_downsample.h5ad", dest = "h5seurat", overwrite = TRUE)
sce.integrated <- LoadH5Seurat("J:/scanpy_out/final_ann/monocle3/FLC_raw_downsample.h5seurat")
sce.integrated <- subset(sce.integrated, subset = celltype%in%c('Telocyte_VSTM2A','Telocyte_NPY','FLC_ACTA2','FLC_SCARA5','FLC_CCL19','FLC_GREM1','FLC_APOD','FLC_CLDN1','FLC_KCNN3'))
DefaultAssay(sce.integrated) <- 'RNA'
data <- GetAssayData(sce.integrated, assay = 'RNA', slot = 'counts')
cell_metadata <- sce.integrated@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
#cds <- align_cds(cds, alignment_group = "ID")
cds <- reduce_dimension(cds, reduction_method = 'UMAP',
                        umap.min_dist = 0.4,
                        umap.n_neighbors = 70)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "partition",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype", label_cell_groups=T,cell_size = 1) +
  scale_color_manual(values=c("#ADCBDA","#3274A0","#B3D395","#41923B","#EFA6A5","#CA3335","#EABC80","#DE7F20","#C7B5D1"))
saveRDS(cds, "J:/scanpy_out/final_ann/monocle3/FLC_traje3.RDS")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- read.csv('J:/scanpy_out/final_ann/monocle3/FLC_UMAP.csv', row.names = 1)
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
#cds <- cluster_cells(cds)
cds <- cluster_cells(cds,k = 20,num_iter = 2,partition_qval = 0.05,resolution = 0.002)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) 
cds <- learn_graph(cds)
cds <- order_cells(cds)
#cds_sub <- choose_graph_segments(cds)
cds_subset <- choose_cells(cds)
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset_pr_test_res[subset_pr_test_res$q_value<0.05 & subset_pr_test_res$morans_I > 0.1,])
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
gene_module_df <- as.data.frame(gene_module_df)
write.csv(gene_module_df, 'J:/scanpy_out/final_ann/pseudotime_CLDN1.csv')
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), 
                                cell_group=colData(cds_subset)$celltype)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
agg_mat <- agg_mat[,c('FLC_GREM1','FLC_APOD','FLC_CLDN1')]
pheatmap(agg_mat,
         scale="column", clustering_method="ward.D2", border = F,
         color = colorRampPalette(c("#2A66A3",'#7AAECC',"#BED8E6",'#F6F5F2',"#F3C2A9",'#D5725C','#A91F2A'))(50))
agg_mat <- as.data.frame(agg_mat)
pheatmap(agg_mat,
         scale="row", clustering_method="ward.D2", border = F,
         color = colorRampPalette(c("#2A66A3",'#7AAECC',"#BED8E6",'#F6F5F2',"#F3C2A9",'#D5725C','#A91F2A'))(50))
genes <- row.names(subset(subset_pr_test_res, q_value < 0.05 & morans_I > 0.1))
pt.matrix <- exprs(cds_subset)[match(genes,rownames(rowData(cds_subset))),order(pseudotime(cds_subset))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "RdBu"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm = draw(htkm)


