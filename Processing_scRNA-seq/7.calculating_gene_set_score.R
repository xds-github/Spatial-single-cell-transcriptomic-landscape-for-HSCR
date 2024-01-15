# convert the h5
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
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
pair11 <- c('#C294B9','#E2D5E6','#990066','#6A371B','#B17D30','#2E796D','#B4B5B5','#C37284','#333366','#D6A128')
# Dealing vessel endothelium
count_matrix <- read.csv("J:/scanpy_out/final_ann/Peri_Endo_out/Endo_peri_raw_matrix_230921.csv",check.names = F)
rownames(count_matrix) <- count_matrix[,1]
count_matrix[,1] <- NULL
count_matrix <- t(count_matrix)
meta_data <- read.csv("J:/scanpy_out/final_ann/Peri_Endo_out/Endo_peri_meta_data.csv",check.names = F)
rownames(meta_data) <- meta_data[,1]
meta_data[,1] <- NULL
sce <- CreateSeuratObject(counts = count_matrix, meta.data = meta_data)
sce <- subset(sce, subset = celltype%in%c("EC_Arterial","EC_Venous"))
sce <- SCTransform(sce, verbose = FALSE)
sce <- RunPCA(sce, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:30, verbose = FALSE)
Idents(sce) <- sce$group
new.cluster.ids <- c('CT', 'HG', 'HA','HG','HA')
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce$group2 <- Idents(sce)
# Zinc score
GO_DATA <- readRDS("D:/Linux_resource/GO_DATA.RDS")
gene_list <- GO_DATA$PATHID2EXTID$'GO:0071294' #Zinc
gene_list <- list(gene_list)
sce <- AddModuleScore(sce,features = gene_list,name = "zinc_score")
VlnPlot(sce, features = c('zinc_score1'),pt.size = 0,cols = c('#759B9C','#E7D3B8','#EA937F')) + geom_boxplot(width = 0.1, fill = 'white')
temp <- sce@meta.data[,c('group2','zinc_score1')]
wilcox.test(temp[temp$group2=='CT',2], temp[temp$group2=='HG',2])
wilcox.test(temp[temp$group2=='CT',2], temp[temp$group2=='HA',2])
wilcox.test(temp[temp$group2=='HG',2], temp[temp$group2=='HA',2])
saveRDS(sce, "J:/scanpy_out/final_ann/Peri_Endo_out/Endo_score.RDS")

# Dealing Pericyte
count_matrix <- read.csv("J:/scanpy_out/final_ann/Peri_Endo_out/Endo_peri_raw_matrix_230921.csv",check.names = F)
rownames(count_matrix) <- count_matrix[,1]
count_matrix[,1] <- NULL
count_matrix <- t(count_matrix)
meta_data <- read.csv("J:/scanpy_out/final_ann/Peri_Endo_out/Endo_peri_meta_data.csv",check.names = F)
rownames(meta_data) <- meta_data[,1]
meta_data[,1] <- NULL
#meta_data <- meta_data[,c("ID","group","batch","clusters","celltype",'group2')]
sce <- CreateSeuratObject(counts = count_matrix, meta.data = meta_data)
sce <- subset(sce, subset = celltype%in%c("Pericyte_RGS5","Pericyte_MYH11","Pericyte_CCL2"))
sce <- SCTransform(sce, verbose = FALSE)
sce <- RunPCA(sce, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:30, verbose = FALSE)
Idents(sce) <- sce$group
new.cluster.ids <- c('CT', 'HG', 'HA','HG','HA')
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce$group2 <- Idents(sce)
sce$celltype <- factor_order_change(c('Pericyte_MYH11','Pericyte_CCL2','Pericyte_RGS5'),sce$celltype)
Idents(sce) <- sce$celltype
# Response to ischemia
gene_list <- GO_DATA$PATHID2EXTID$'GO:0002931' #Response to ischemia
gene_list <- list(gene_list)
sce <- AddModuleScore(sce,features = gene_list,name = "ischemia_score")
VlnPlot(sce, features = c('ischemia_score1'),pt.size = 0,cols = c('#B2D2CE','#D58889','#7594A7')) + geom_boxplot(width = 0.1, fill = 'white')
temp <- sce@meta.data[,c('celltype','ischemia_score1')]
wilcox.test(temp[temp$celltype=='Pericyte_RGS5',2], temp[temp$celltype=='Pericyte_MYH11',2])
wilcox.test(temp[temp$celltype=='Pericyte_RGS5',2], temp[temp$celltype=='Pericyte_CCL2',2])
wilcox.test(temp[temp$celltype=='Pericyte_MYH11',2], temp[temp$celltype=='Pericyte_CCL2',2])
saveRDS(sce, "J:/scanpy_out/final_ann/Peri_Endo_out/Pericyte_score.RDS")
