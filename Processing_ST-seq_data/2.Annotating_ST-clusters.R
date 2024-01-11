#library(SPATA2)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
pair11 <- c('#E64B35','#4DBBD5','#00A087','#3C5488','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
temp <- read.csv("F://data_for_upload/ST-seq/ST-seq_chip_info/ST-seq_meta_info.txt", sep = '\t', header = F)
temp <- temp$V1
for (i in temp) {
  sce2 <- Seurat::Load10X_Spatial(paste0(i,"/outs/"),filter.matrix = F)
  sce2$ID <- i
  sce2 <- subset(sce2, subset = nFeature_Spatial>400)
  sce2 <- SCTransform(sce2, assay = "Spatial", verbose = FALSE)
  sce2 <- RunPCA(sce2, assay = "SCT", verbose = FALSE, npcs = 60)
  sce2 <- FindNeighbors(sce2, reduction = "pca", dims = 1:30)
  saveRDS(sce2, paste0(i,'.RDS'))
}
HD30K <- readRDS('201132A_HD30K.RDS')
HD30K$donor <- 'HD30'
HD30X <- readRDS('201132A_HD30X.RDS')
HD30X$donor <- 'HD30'
HD30Y <- readRDS('201132A_HD30Y.RDS')
HD30Y$donor <- 'HD30'
HD32K <- readRDS('201132A_HD32K.RDS')
HD32K$donor <- 'HD32'
HD32X <- readRDS('201132B_HD32X.RDS')
HD32X$donor <- 'HD32'
HD33K <- readRDS('201132B_HD33K.RDS')
HD33K$donor <- 'HD33'
HD33Y <- readRDS('201132B_HD33Y.RDS') 
HD33Y$donor <- 'HD33'
HD35K <- readRDS('201132C_HD35K.RDS')
HD35K$donor <- 'HD35'
HD35X <- readRDS('201132C_HD35X.RDS')
HD35X$donor <- 'HD35'
HD36K <- readRDS('201132C_HD36K.RDS')
HD36K$donor <- 'HD36'
HD36X <- readRDS('201132C_HD36X.RDS')
HD36X$donor <- 'HD36'
HD36Y <- readRDS('201132C_HD36Y.RDS')
HD36Y$donor <- 'HD36'
HD37K <- readRDS('201132C_HD37K.RDS')
HD37K$donor <- 'HD37'
HD37X <- readRDS('201132C_HD37X.RDS')
HD37X$donor <- 'HD37'
HD37Y <- readRDS('201132C_HD37Y.RDS')
HD37Y$donor <- 'HD37'
HD39K <- readRDS('201132C_HD39K.RDS')
HD39K$donor <- 'HD39'
HD39X <- readRDS('201132C_HD39X.RDS')
HD39X$donor <- 'HD39'
HD39Y <- readRDS('201132C_HD39Y.RDS')
HD39Y$donor <- 'HD39'
CT10 <- readRDS('201132D_CT10.RDS')
CT10$donor <- 'CT10'
CT11 <- readRDS('201132D_CT11.RDS')
CT11$donor <- 'CT11'
CT12 <- readRDS('201132D_CT12.RDS')
CT12$donor <- 'CT12'
HD32Y <- readRDS('201132D_HD32Y.RDS')
HD32Y$donor <- 'HD32'
CT14 <- readRDS('201132E_CT14.RDS')
CT14$donor <- 'CT14'
CT17 <- readRDS('201132E_CT17.RDS')
CT17$donor <- 'CT17'
ST_merge <- merge(HD30K,c(HD30X,HD30Y,HD32K,HD32X,HD33K,HD33Y,HD35K,HD35X,HD36K,HD36X,HD36Y,HD37K,HD37X,HD37Y,HD39K,HD39X,HD39Y,CT10,CT11,CT12,HD32Y,CT14,CT17))
DefaultAssay(ST_merge) <- "SCT"
VariableFeatures(ST_merge) <- c(VariableFeatures(HD30K), VariableFeatures(HD30X), VariableFeatures(HD30Y), 
                                VariableFeatures(HD32K), VariableFeatures(HD32X), VariableFeatures(HD33K), 
                                VariableFeatures(HD33Y), VariableFeatures(HD35K), VariableFeatures(HD35X), 
                                VariableFeatures(HD36K), VariableFeatures(HD36X), VariableFeatures(HD36Y), 
                                VariableFeatures(HD37K), VariableFeatures(HD37X), VariableFeatures(HD37Y), 
                                VariableFeatures(HD39K), VariableFeatures(HD39X), VariableFeatures(HD39Y), 
                                VariableFeatures(CT10), VariableFeatures(CT11), VariableFeatures(CT12), 
                                VariableFeatures(HD32Y), VariableFeatures(CT14), VariableFeatures(CT17))
rm(HD30K,HD30X,HD30Y,HD32K,HD32X,HD33K,HD33Y,HD35K,HD35X,HD36K,HD36X,HD36Y,HD37K,HD37X,HD37Y,HD39K,HD39X,HD39Y,CT10,CT11,CT12,HD32Y,CT14,CT17)
ST_merge <- RunPCA(ST_merge, verbose = FALSE)
ST_merge <- FindNeighbors(ST_merge, dims = 1:30)
ST_merge <- FindClusters(ST_merge, verbose = FALSE)
ST_merge <- RunUMAP(ST_merge, dims = 1:30)
DimPlot(ST_merge, label = T)
DimPlot(ST_merge, group.by = 'celltype', label = T)
ST_merge$SCT_snn_res.0.8 <- ST_merge$SCT_snn_res.1.2 <- ST_merge$SCT_snn_res.1 <- ST_merge$SCT_snn_res.1.5 <- ST_merge$SCT_snn_res.1.7 <- ST_merge$SCT_snn_res.2 <- ST_merge$SCT_snn_res.1.3 <- ST_merge$SCT_snn_res.1.6 <- NULL
ST_merge$fun_cluster <- ST_merge$celltype
#---------------------------------------- Annotate different functional region ---------------------------------------------
Idents(ST_merge) <- ST_merge$seurat_clusters
ST_merge <- set_celltype(ST_merge, new.cluster.ids = c('C6','C5','C9','C11','C2','C8','C2','C1','C8','C8','C3','C12',
                                                       'C8','C9','C7','C4','C10','C9','C8','C9','C13','C9','C9','C8',
                                                       'C9','C3','C9','C3','C9','C8','C8'))
ST_merge$celltype <- factor_order_change(c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13'), ST_merge$celltype)
ST_merge$fun_cluster2 <- ST_merge$celltype
Idents(ST_merge) <- ST_merge$fun_cluster2
DimPlot(ST_merge, label = T, cols = pair11, pt.size = 1)

Idents(ST_merge) <- ST_merge$ID
ST_merge <- set_celltype(ST_merge, new.cluster.ids = c("HG","HA","HA","HG","HA","HG","HG","HG",
                                                       "HA","HG","HA","HA","HG","HA","HA","HG",
                                                       "HA","HG","CT","CT","CT","HG","CT","CT"))
ST_merge$celltype <- factor_order_change(c('CT','HG','HA'), ST_merge$celltype)
ST_merge$group <- ST_merge$celltype
ST_merge$celltype <- NULL
saveRDS(ST_merge, 'J:/ST-analysis/total_merge2.RDS')
#------------------------------------------------ Plot marker genes ------------------------------------------------------------
avexpress_RNA <- AverageExpression(ST_merge)
gene_list <- c('EPCAM','BEST4','ZG16','ACHE','SLC26A3',
               'MUC2','MUC1','IGHA2','IGHA1','ADAMDEC1','OLFM4','CXCL14',
               'GREM1','PI15','AEBP1',
               'PTPRC','CD3E','CD19','CD69','AICDA','IGHG4','CCL19','LAMP3','CCR7',
               'C3','CFD','COL1A1','COL1A2','FBLN1',
               'MALAT1','MYH11','ACTG2',
               'ACTA2','JUN','ZFP36','SLC30A1','TNFRSF12A',
               'FOSL1','CCL2','CXCL8','THBS1','ICAM1',
               'NOTCH3','PECAM1','HBB',
               'VIP','GAL','NOS1','MAP2','ALDH1A1','BCHE','S100B',
               'GPM6B','APOD','CLDN1')
avexpress_RNA <- avexpress_RNA$SCT
avexpress_RNA <- avexpress_RNA[gene_list,]
library(pheatmap)
avexpress_RNA1 <- t(scale(t(avexpress_RNA), scale=T, center = T))
pheatmap(avexpress_RNA1, border = F, show_rownames = T, cluster_cols = F, color = colorRampPalette(c("#2A66A3",'#7AAECC',"#BED8E6",'#F6F5F2',"#F3C2A9",'#D5725C','#A91F2A'))(50), cluster_rows = F)
#----------------------------------------------- Plot distribution --------------------------------------------------------
library(reshape2)
temp <- table(ST_merge$ID,ST_merge$fun_cluster2)
#write.csv(temp,'ST-cluster_distr_processed.csv')
celltype_count <- read.csv('ST-cluster_distr_processed.csv')
celltype_count$donor <- celltype_count$X
celltype_count$X <- NULL
celltype_count <- melt(celltype_count, id.vars = c("donor"))
celltype_count$donor <- factor_order_change(c('201132E_CT14','201132D_CT11','201132E_CT17','201132D_CT10','201132D_CT12',
                                              '201132A_HD30K','201132C_HD35K','201132B_HD33G-2','201132C_HD36K','201132C_HD37K',
                                              '201132C_HD39K','201132A_HD32K','201132D_HD32G-2','201132B_HD33Y',
                                              '201132C_HD39Y','201132C_HD36Y','201132A_HD30Y','201132C_HD37Y','201132C_HD35X',
                                              '201132B_HD32X','201132C_HD39X','201132C_HD37X','201132A_HD30X','201132C_HD36X'),celltype_count$donor)
ggplot(celltype_count, aes(donor, value, fill = variable)) + theme_classic() +
  geom_bar(stat='identity',position='fill', colour = "white") + scale_fill_manual(values = pair11) + RotatedAxis()
#----------------------------------------------- Obtain DEGs --------------------------------------------------------
temp <- FindAllMarkers(ST_merge)
write.csv(temp, 'ST-cluster_total_marker.csv')

temp <- FindMarkers(ST_merge,ident.1 = 'C13', ident.2 = 'C12', logfc.threshold = 0)
write.csv(temp, 'J:/ST-analysis/C13_v_C12.csv')

temp <- FindMarkers(ST_merge,ident.1 = 'C10', ident.2 = c('C8','C9'), logfc.threshold = 0)
write.csv(temp, 'J:/ST-analysis/C10_v_restCMusc.csv')

#------------------------------------------------ Exporting the meta data ------------------------------------------------------------
temp <- ST_merge@meta.data
barcode <- substring(rownames(temp), 1,18)
temp$trans_barcode <- paste(temp$ID, barcode, sep = '_')
write.csv(temp, 'J:/ST-analysis/ST_metadata.csv')
#------------------------------------------------ Visualizing the representative data ------------------------------------------------------------
pair11 <- c('#E64B35','#4DBBD5','#00A087','#3C5488','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.17', pt.size.factor = 1.2, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = F)

pair11 <- c('#E64B35','#4DBBD5','#00A087','#3C5488','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#7E6148','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.16', pt.size.factor = 1.2, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = F)


pair11 <- c('#E64B35','#4DBBD5','#00A087','#3C5488','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.20', pt.size.factor = 1.5, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = F)

pair11 <- c('#E64B35','#4DBBD5','#00A087','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.15', pt.size.factor = 1.2, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = F)

pair11 <- c('#E64B35','#4DBBD5','#00A087','#3C5488','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.16', pt.size.factor = 1.2, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = F)

pair11 <- c('#E64B35','#4DBBD5','#00A087','#3C5488','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.12', pt.size.factor = 1.2, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = T)

pair11 <- c('#E64B35','#4DBBD5','#00A087','#8491B4','#DCDDDD','#F39B7F','#91D1C2','#0099B4','#FEB500','#7E6148','#EF7000','#643C90')
SpatialDimPlot(ST_merge, images = 'slice1.15', pt.size.factor = 1.5, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11, crop = T)
SpatialPlot(ST_merge, images = 'slice1.15', pt.size.factor = 1.5, features = 'HSPA1A') +
  scale_fill_gradient2(low = "#2165AB", mid = "white", high = "#E64B35", midpoint = 2)
SpatialPlot(ST_merge, images = 'slice1.15', pt.size.factor = 1.5, features = 'CCL2') +
  scale_fill_gradient2(low = "#2165AB", mid = "white", high = "#E64B35", midpoint = 1.7)
SpatialPlot(ST_merge, images = 'slice1.15', pt.size.factor = 1.5, features = 'SLC30A1') +
  scale_fill_gradient2(low = "#2165AB", mid = "white", high = "#E64B35", midpoint = 1.7)
SpatialPlot(ST_merge, images = 'slice1.15', pt.size.factor = 1.5, features = 'FOSL1') +
  scale_fill_gradient2(low = "#2165AB", mid = "white", high = "#E64B35", midpoint = 0.8)
SpatialFeaturePlot(ST_merge, images = 'slice1.16', pt.size.factor = 1.2, features = 'S100B', alpha = c(0.5, 1))
SpatialPlot(ST_merge, images = 'slice1.16', pt.size.factor = 1.2, features = 'APOD', alpha = c(0.1, 1))+ scale_fill_gradient(low="#3295F7", high="#950C00")
SpatialPlot(ST_merge, images = 'slice1.16', pt.size.factor = 1.5, features = 'APOD', alpha = c(0.1, 1))+ scale_fill_distiller(palette = "RdBu")
Idents(ST_merge) <- ST_merge$fun_cluster2
SpatialDimPlot(ST_merge, images = 'slice1.2', pt.size.factor = 1.5, label = T,stroke = 0, group.by = 'fun_cluster2', cols = pair11)
SpatialDimPlot(ST_merge, images = 'slice1.2', pt.size.factor = 1.5, label = T,stroke = 0, group.by = 'seurat_clusters')
VlnPlot(ST_merge, features = 'nFeature_Spatial', pt.size = 0, cols = pair11) + geom_boxplot(width = 0.1, fill = 'white')
DotPlot(ST_merge, features = c('NLRP1','NLRP3','NLRP12','NLRP6','NLRP7','NLRC4','NLRP4','AIM2','IFI16','PYCARD','CASP1','GSDMD','SCAF11','CASP4','CASP5','CASP2')) +RotatedAxis()

#------------------------------------------------ Calculating gene set score  ------------------------------------------------------------
st_temp <- subset(ST_merge, subset = fun_cluster2%in%c('C8','C9','C10'))
GO_DATA <- readRDS("D:/Linux_resource/GO_DATA.RDS")
gene_list <- GO_DATA$PATHID2EXTID$'GO:0002931'
gene_list <- list(gene_list)
st_temp <- AddModuleScore(st_temp,features = gene_list,name = "is_score")
Idents(st_temp) <- st_temp$fun_cluster2
VlnPlot(st_temp, features = c('is_score1'),cols = c('#91D1C2','#0099B4','#FEB500'),pt.size = 0) + geom_boxplot(width = 0.1, fill = 'white')
temp <- st_temp@meta.data[,c('fun_cluster2','is_score1')]
wilcox.test(temp[temp$fun_cluster2=='C10','is_score1'],temp[temp$fun_cluster2=='C9','is_score1'])
wilcox.test(temp[temp$fun_cluster2=='C10','is_score1'],temp[temp$fun_cluster2=='C8','is_score1'])
saveRDS(st_temp, "J:/ST-analysis/Musclar_ST_clusters_score.RDS")

#------------------------------------------------ Extract raw matrix for cibersortx --------------------------------------------------------
library(Matrix)
ST_merge <- subset(ST_merge, subset = fun_cluster2!='C6')
ST_matrix <- ST_merge@assays$Spatial@counts
ST_matrix <- as.matrix(ST_matrix)
ST_matrix <- t(ST_matrix)
rownames(ST_matrix) <- ST_merge@meta.data$fun_cluster2
ST_matrix <- t(ST_matrix)
ST_matrix <- ST_matrix[rowSums(ST_matrix) > 0,]
m_colname <-colnames(ST_matrix)
m_colname <- c('GeneSymbol',m_colname)
ST_matrix <- data.frame(ST_matrix)
ST_matrix <- ST_matrix %>%  mutate(GeneSymbol = rownames(ST_matrix)) %>% select(GeneSymbol, everything())
colnames(ST_matrix) <- m_colname
write.table(ST_matrix, "ST_raw_matrix.txt", row.names = F, quote = F)
