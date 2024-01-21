# This script was modified from  Gary Reynolds et al.'s work. Please cite Reynolds, G. et al. Developmental cell programs are co-opted in inflammatory skin disease. Science 371 (2021). 
library(MASS)
library(reshape2)
count_matrix <- read.csv('/share/home/xudeshu/scanpy_dic/HSCR/final_anan/meta_file/total_meta_data2.csv')
count_matrix <- melt(count_matrix)
cond_1 <- "HA"
cond_2 <- "HG"
cond_3 <- "CT"
pv_cond_2vs1 <- vector()
pv_cond_3vs1 <- vector()
pv_cond_3vs2 <- vector()
# Compute the p-values for increasing labeled fraction (using negative binomial regression)
for(i in 1:length(unique(count_matrix$variable))) {
  cell = levels(count_matrix$variable)[i]
  totals = aggregate(count_matrix$value, by = list(count_matrix$donor), FUN = sum)
  colnames(totals) <- c("donor", "total")
  d = data.frame(subset(count_matrix, variable == cell))
  d <- merge(d, totals, by = "donor")
  
  if(any(d$total == 0)) {
    warning("Removing sample that didn't detect this celltype")
    d = subset(d, total > 0)
  }
  
  { # Copy this block as many times as many condition-pairs you test
    reference <- cond_1 # !
    test <- cond_2 # !
    
    x1 = d[d$group %in% c(reference, test), ]
    nb = MASS::glm.nb(formula = value ~ group + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
    pv_cond_2vs1[i] = anova(nb, test = "LRT")$`Pr(>Chi)`[2] # !
  }
  
  { # Copy this block as many times as many condition-pairs you test
    reference <- cond_1 # !
    test <- cond_3 # !
    
    x1 = d[d$group %in% c(reference, test), ]
    nb = MASS::glm.nb(formula = value ~ group + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
    pv_cond_3vs1[i] = anova(nb, test = "LRT")$`Pr(>Chi)`[2] # !
  }
  
  { # Copy this block as many times as many condition-pairs you test
    reference <- cond_2 # !
    test <- cond_3 # !
    
    x1 = d[d$group %in% c(reference, test), ]
    nb = MASS::glm.nb(formula = value ~ group + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
    pv_cond_3vs2[i] = anova(nb, test = "LRT")$`Pr(>Chi)`[2] # !
  }
}
pv = cbind.data.frame(pv_cond_2vs1, pv_cond_3vs1,pv_cond_3vs2) # !
rownames(pv) = levels(count_matrix$variable)
write.csv(pv, file = "J:/scanpy_out/final_ann/pv_celltount_all.csv")


