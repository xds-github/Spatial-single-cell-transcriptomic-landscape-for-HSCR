library(Seurat)
library(dplyr)
library(patchwork)
library(pheatmap)
factor_order_change <- function(new_order, old_factor){
  new_order <- factor(1:length(new_order),labels = new_order)
  new_factor <- factor(old_factor,levels = levels(new_order))
  return(new_factor)
}
set_celltype <- function(sce, new.cluster.ids = c(1,2)){
  names(new.cluster.ids) <- levels(sce)
  sce <- RenameIdents(sce, new.cluster.ids)
  sce$celltype <- Idents(sce)
  return(sce)
}
sce.integrated <- readRDS("J:/scanpy_out/final_ann/cell_interaction/total_reference_raw2.RDS")
sce.integrated$celltype <- factor_order_change(c('Myocyte','Telocyte_VSTM2A','Telocyte_NPY','FLC_ACTA2','FLC_SCARA5','FLC_CCL19','FLC_GREM1','FLC_APOD','FLC_CLDN1','FLC_KCNN3',
                                                 'Pericyte_RGS5','Pericyte_CCL2','Pericyte_MYH11','EC_Lym','EC_Venous','EC_Arterial',
                                                 'Epi_stem', 'TA','Enterocytes', 'Enterocytes_SLC26A3','Enterocytes_LYZ', 'Enterocytes_BEST4', 'Goblet_RETNLB', 'Goblet_ACHE',  'Tuft', 'Enteroendocrine',
                                                 'Glial','Neuron',
                                                 'GCB_MKI67','GCB','ProB','Brm_IGL','Brm_FCRL4','Brm_CD267','Brm','Plasmablast','Plasma',
                                                 'LTi','Th17','Th2','Treg_like','Tfh_like','CD4T_naive','CD8T_naive','CD8Tcm','CD8Tgd_ZNF683','CD8Tgd_CD39','CD8Trm','CD8Tem','CD8Tem_CX3CR1','T_cycling','gdT','NK',
                                                 'M_IL1B','M_C1Q','M_LYVE1','M_APOE','M_cycling','cDC2','cDC1','DC_LAMP3','pDC','Mast'), sce.integrated$celltype)
Idents(sce.integrated) <- sce.integrated$celltype
sce.integrated <- SCTransform(sce.integrated, verbose = FALSE)
sce.integrated <- RunPCA(sce.integrated, verbose = FALSE)
sce.integrated <- RunUMAP(sce.integrated, dims = 1:30, verbose = FALSE)
ST_merge <- readRDS('ST-total_merge.RDS')
ST_merge <- subset(ST_merge, fun_cluster2%in%c('C1','C2','C3','C4','C5','C7','C8','C9','C10','C11','C12','C13'))
Idents(ST_merge) <- ST_merge$fun_cluster2
anchors <- FindTransferAnchors(reference = sce.integrated, query = ST_merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = sce.integrated$celltype, prediction.assay = TRUE,
                                  weight.reduction = ST_merge[["pca"]], dims = 1:30)
ST_merge[["predictions"]] <- predictions.assay
DefaultAssay(ST_merge) <- "predictions"
#celltype_list
celltype_list <- c('FLC-ACTA2','FLC-SCARA5','FLC-CCL19','FLC-GREM1','FLC-APOD','FLC-CLDN1','FLC-KCNN3',
                   'Pericyte-RGS5','Pericyte-MYH11','Pericyte-CCL2','EC-Lym','EC-Venous','EC-Arterial',
                   'Epi-stem', 'TA','Enterocytes', 'Enterocytes-SLC26A3','Enterocytes-LYZ', 'Enterocytes-BEST4', 'Goblet-RETNLB', 'Goblet-ACHE',
                   'Glial','Neuron',
                   'GCB-MKI67','GCB','ProB','Brm-IGL','Plasmablast','Plasma',
                   'M-C1Q','M-LYVE1','M-APOE','M-cycling','cDC1','DC-LAMP3','pDC','Mast')

DotPlot(ST_merge, features = celltype_list, cols = 'RdBu', scale.by = 'radius') + RotatedAxis()

# Abundance analysis
ST_temp <- subset(ST_merge, subset = fun_cluster2%in%c('C8','C9','C10'))
VlnPlot(ST_temp, features = c('Pericyte-MYH11'), pt.size = 0, cols = c('#91D1C2','#0099B4','#FEB500')) + geom_boxplot(width = 0.1, fill = 'white')
VlnPlot(ST_temp, features = c('Pericyte-CCL2'), pt.size = 0, cols = c('#91D1C2','#0099B4','#FEB500')) + geom_boxplot(width = 0.1, fill = 'white')
VlnPlot(ST_temp, features = c('FLC-ACTA2'), pt.size = 0, cols = c('#91D1C2','#0099B4','#FEB500')) + geom_boxplot(width = 0.1, fill = 'white')
VlnPlot(ST_temp, features = c('FLC-KCNN3'), pt.size = 0, cols = c('#91D1C2','#0099B4','#FEB500')) + geom_boxplot(width = 0.1, fill = 'white')

ST_temp <- subset(ST_merge, subset = fun_cluster2%in%c('C12','C13'))
VlnPlot(ST_temp, features = c('Neuron'), pt.size = 0, cols = c('#EF7000','#643C90')) + geom_boxplot(width = 0.1, fill = 'white')
VlnPlot(ST_temp, features = c('Glial'), pt.size = 0, cols = c('#EF7000','#643C90')) + geom_boxplot(width = 0.1, fill = 'white')
VlnPlot(ST_temp, features = c('FLC-APOD'), pt.size = 0, cols = c('#EF7000','#643C90')) + geom_boxplot(width = 0.1, fill = 'white')
VlnPlot(ST_temp, features = c('FLC-CLDN1'), pt.size = 0, cols = c('#EF7000','#643C90')) + geom_boxplot(width = 0.1, fill = 'white')
