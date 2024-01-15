library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
library(edgeR)
#----------------------------------------- Generating ST_reference -----------------------------------------------------------
# Down sample  each cluster to 500
ST_merge <- readRDS('J:/ST-analysis/ST-total_merge.RDS')
temp <- ST_merge@meta.data
temp2 <- c()
for (i in c('C1','C2','C3','C4','C5','C7','C8','C9','C10','C11','C12','C13')) {
  temp1 <- temp[temp$fun_cluster2==i,]
  temp1 <- rownames(temp1)
  temp1 <- sample(temp1, 500)
  temp2 <- c(temp2, temp1)
}
temp <- temp[temp2,]
ST_matrix <- ST_merge@assays$Spatial@counts
ST_matrix <- as.matrix(ST_matrix)
ST_matrix <- ST_matrix[rowSums(ST_matrix) > 0,]
ST_matrix <- ST_matrix[,rownames(temp)]
colnames(ST_matrix) <- temp$fun_cluster2
ST_matrix <- data.frame(ST_matrix)
ST_matrix <- ST_matrix %>%  mutate(GeneSymbol = rownames(ST_matrix)) %>% select(GeneSymbol, everything())
#colnames(ST_matrix) <- c('GeneSymbol',as.character(temp$fun_cluster2))
#ST_matrix <- as.matrix(ST_matrix)
write.table(ST_matrix,"J:/ST-analysis/ST_raw_matrix.txt", row.names = F, quote = F, sep = '\t')
