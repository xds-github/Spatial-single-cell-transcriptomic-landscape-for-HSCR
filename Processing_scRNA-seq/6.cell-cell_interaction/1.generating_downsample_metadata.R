library(Seurat)
factor_order_change <- function(new_order, old_factor){
  new_order <- factor(1:length(new_order),labels = new_order)
  new_factor <- factor(old_factor,levels = levels(new_order))
  return(new_factor)
}
target_mata <- read.csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/total_meta_data2.csv')
celltype <- unique(target_mata$celltype)
#down sample the metadata
filter_meta <- data.frame()
for (i in celltype) {
  temp <- target_mata[target_mata$celltype==i,]
  if (dim(temp)[1]>500) {
    temp1 <- sample(temp$X, 500)
    temp2 <- temp[temp$X%in%temp1,]
    filter_meta <- rbind(filter_meta, temp2)
  } else {
    filter_meta <- rbind(filter_meta, temp)
  }
}
filter_meta$X <- NULL
write.csv(filter_meta, '/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/subset_total_500_meta.csv',quote = F, row.names = F)
