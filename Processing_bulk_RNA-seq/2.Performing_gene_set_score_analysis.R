library(DESeq2)
library(apeglm)
library(IHW)
library(org.Hs.eg.db)
library(biomaRt)
library(curl)
library(ggplot2)
library(Seurat)
library(dplyr)
library(edgeR)
set_celltype <- function(sce, new.cluster.ids = c(1,2)){
  names(new.cluster.ids) <- levels(sce)
  sce <- RenameIdents(sce, new.cluster.ids)
  sce$celltype <- Idents(sce)
  return(sce)
}
factor_order_change <- function(new_order, old_factor){
  new_order <- factor(1:length(new_order),labels = new_order)
  new_factor <- factor(old_factor,levels = levels(new_order))
  return(new_factor)
}
# ----------------------- Read in matrix -----------------------------------
meta_data <- read.csv('Bulk_RNA-seq_info.csv',check.names = F)
count_matrix <- read.table(paste0('E:/megacolon/Bulk-RNAseq/star_out_220429/','20109903X','_star_out.txt'))
count_matrix$V2 <- count_matrix$V3 <- NULL
for (i in sample_N) {
  temp <- read.table(paste0('E:/megacolon/Bulk-RNAseq/star_out_220429/',i,'_star_out.txt'))
  temp <- plyr::rename(temp, c('V4'=i))
  temp$V2 <- temp$V3 <- NULL
  count_matrix <- merge(count_matrix, temp, by = 'V1')
}
rownames(count_matrix) <- count_matrix$V1
count_matrix$V1 <- count_matrix$V4 <- NULL
# ----------------------- transform ID-----------------------------------
refer_gene_list <- read.table('D:/Linux_resource/Ensembl_human_mouse_ID_change.txt', sep = '\t', header = T,check.names = F)
refer_gene_list <- refer_gene_list[refer_gene_list$`Gene stable ID`%in%rownames(count_matrix),]
refer_gene_list <- refer_gene_list[!duplicated(refer_gene_list$`Gene stable ID`),]
dup_list <- unique(refer_gene_list$`Gene name`[duplicated(refer_gene_list$`Gene name`)])
dup_list <- refer_gene_list[refer_gene_list$`Gene name`%in%dup_list,'Gene stable ID']
refer_gene_list$`Gene_name2` <- paste(refer_gene_list$`Gene stable ID`,refer_gene_list$`Gene name` ,sep = '_')
for (i in dup_list) {
  refer_gene_list[refer_gene_list$`Gene stable ID`==i,'Gene name'] <-  refer_gene_list[refer_gene_list$`Gene stable ID`==i,'Gene_name2']
}
refer_gene_list$`Gene_name2` <- NULL
count_matrix <- count_matrix[refer_gene_list$`Gene stable ID`,]
rownames(count_matrix) <- refer_gene_list$`Gene name`
#saveRDS(count_matrix, 'E:/megacolon/Bulk-RNAseq/star_out_220429/total_raw_count_matrix_8.RDS')
# ----------------------- Generating sce object ----------------------------------
sce <- CreateSeuratObject(counts = count_matrix)
sce$group <- meta_data$group
sce$segment <- meta_data$segment
sce$batch <- meta_data$batch
sce$ID <- meta_data$ID
sce <- PercentageFeatureSet(sce, pattern = "^MT-", col.name = "percent.mt")
sce <- SCTransform(sce)
sce <- RunPCA(sce)
Idents(sce) <- sce$segment
# caculate apoptosis score
GO_DATA <- readRDS("D:/Linux_resource/GO_DATA.RDS")
gene_list <- GO_DATA$PATHID2EXTID$'GO:0043525'
write.csv(gene_list, "J:/scanpy_out/final_ann/GO_neural_apop_list.csv")
gene_list <- list(gene_list)
sce <- AddModuleScore(sce,features = gene_list,name = "apoptosis_score")
VlnPlot(sce, features = c('apoptosis_score1'))#positive regulation of neuron apoptotic process
temp <- sce@meta.data[,c('segment','apoptosis_score1')]
wilcox.test(temp[temp$segment=='HA','apoptosis_score1'],temp[temp$segment=='CT','apoptosis_score1'])
wilcox.test(temp[temp$segment=='HG','apoptosis_score1'],temp[temp$segment=='CT','apoptosis_score1'])
wilcox.test(temp[temp$segment=='HG','apoptosis_score1'],temp[temp$segment=='HA','apoptosis_score1'])
temp %>% mutate(segment = factor(segment, levels = c('CT','HG','HA'))) %>%
  ggplot(aes(segment,apoptosis_score1)) + geom_violin(scale = 'width', aes(fill = segment)) +
  geom_boxplot(width = 0.1, fill = 'white') +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA)) +
  scale_fill_manual(values=c('#759B9C','#E7D3B8','#EA937F'))

temp %>% mutate(segment = factor(segment, levels = c('CT','HG','HA'))) %>%
  ggplot(aes(segment,apoptosis_score1, fill = segment)) + geom_boxplot(width = 0.5, fill = 'white') +
  geom_jitter(aes(fill  = segment, color = segment),width = 0.3, size = 2)+
  scale_color_manual(values=c('#759B9C','#E7D3B8','#EA937F')) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA))
#saveRDS(sce, "E:/megacolon/Bulk-RNAseq/rmats-analysis/total_merge8.RDS")
# Acquire CPM matrix
count_matrix <- readRDS('E:/megacolon/Bulk-RNAseq/star_out_220429/total_raw_count_matrix_7.RDS')
count_matrix <- count_matrix[rowSums(count_matrix)>0,]
CPM_matrix <- cpm(count_matrix)
gene <- rownames(CPM_matrix)
CPM_matrix <- cbind(gene, CPM_matrix)
saveRDS(CPM_matrix, 'E:/megacolon/Bulk-RNAseq/star_out_220429/CPM_count_matrix_7.RDS')
write.table(CPM_matrix, 'E:/megacolon/Bulk-RNAseq/star_out_220429/CPM_count_matrix_7.txt', sep = '\t', quote = F, row.names = F)
