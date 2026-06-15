#BiocManager::install("GSVA")
library(GSVA) 
library(limma)
library(msigdbr)
library(dplyr)
library(Seurat)
library(BiocParallel)
library(ggplot2)

GO_DATA <- readRDS("D:/Linux_resource/GO_DATA.RDS")
#GO_LIST <- as.list(GO_DATA)
sce <- readRDS("E:/HSCR/STM_subm/Object_of_bulk_RNA-seq.RDS")
meta_info <- sce@meta.data
dat <- as.matrix(sce@assays$RNA$counts)
#dat <- as.matrix(log2(tpm+1))
#gsva_mat <- gsva(expr=dat, 
#                 gset.idx.list=GO_DATA$PATHID2EXTID, 
#                 kcdf="Gaussian",#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
#                 verbose=T,
#                 parallel.sz = parallel::detectCores())
gsva_mat <- gsva(expr=dat, 
                 gset.idx.list=GO_DATA$PATHID2EXTID, 
                 kcdf="Poisson",#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T,
                 parallel.sz = parallel::detectCores())
gsva_mat <- gsva_mat[,rownames(meta_info)]
#write.csv(gsva_mat,"E:/HSCR/STM_subm/bulk_gsva_go_matrix.csv")
# 可视化ischemia的结果
#gsva_mat <- read.csv("E:/HSCR/STM_subm/bulk_gsva_go_matrix.csv", row.names = 1)

is_core <- gsva_mat['GO:0002931',]
#is_core <- is_core[,paste0('X',rownames(meta_info))]
is_core <- is_core[rownames(meta_info)]
meta_info$is_core1 <- is_core
meta_info$GO_0002931 <- is_core
meta_info %>% mutate(segment = factor(segment, levels = c('CT','HG','HA'))) %>%
  ggplot(aes(segment,is_core1, fill = segment)) + geom_boxplot(width = 0.5, fill = 'white') +
  geom_jitter(aes(fill  = segment, color = segment),width = 0.3, size = 2)+
  scale_color_manual(values=c('#759B9C','#E7D3B8','#EA937F')) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA))

# 可视化血小板聚集的结果
is_core <- gsva_mat['GO:0047179',]
is_core <- is_core[rownames(meta_info)]
meta_info$is_core1 <- is_core
meta_info$GO_0047179 <- is_core
meta_info %>% mutate(segment = factor(segment, levels = c('CT','HG','HA'))) %>%
  ggplot(aes(segment,is_core1, fill = segment)) + geom_boxplot(width = 0.5, fill = 'white') +
  geom_jitter(aes(fill  = segment, color = segment),width = 0.3, size = 2)+
  scale_color_manual(values=c('#759B9C','#E7D3B8','#EA937F')) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = 'black',fill = NA))
meta_info$sample_ID <- rownames(meta_info)
temp <- meta_info[,c('sample_ID','segment','GO_0002931','GO_0047179')]
write.csv(temp, "E:/HSCR/STM_subm/gsva_score_in_bulk.csv", row.names = F)
# 比较HG和CT
submeta <- meta_info[meta_info$segment%in%c('CT','HG'),]
subgsva <- gsva_mat[,rownames(submeta)]
Group <- factor(submeta$segment,levels = c("HG","CT"))
design <- model.matrix(~Group)
fit <- lmFit(subgsva,design)
fit <- eBayes(fit)
tempOutput <- topTable(fit, n=Inf)
tempOutput['GO:0047179',]
tempOutput['GO:0002931',]
# 比较HA和CT
submeta <- meta_info[meta_info$segment%in%c('CT','HA'),]
subgsva <- gsva_mat[,rownames(submeta)]
Group <- factor(submeta$segment,levels = c("HA","CT"))
design <- model.matrix(~Group)
fit <- lmFit(subgsva,design)
fit <- eBayes(fit)
tempOutput <- topTable(fit, n=Inf)
tempOutput['GO:0047179',]
tempOutput['GO:0002931',]
# 比较HA和HG
submeta <- meta_info[meta_info$segment%in%c('HG','HA'),]
subgsva <- gsva_mat[,rownames(submeta)]
Group <- factor(submeta$segment,levels = c("HA","HG"))
design <- model.matrix(~Group)
fit <- lmFit(subgsva,design)
fit <- eBayes(fit)
tempOutput <- topTable(fit, n=Inf)
tempOutput['GO:0047179',]
tempOutput['GO:0002931',]
