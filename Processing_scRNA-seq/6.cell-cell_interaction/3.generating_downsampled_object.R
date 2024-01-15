#--------------------------------------------------------------------------------------- Generating a downsampled object ------------------------------------------------------------------------------
library(rhdf5)
library(Matrix)
# import scRNA-seq data
mydata <- h5read("J:/scanpy_out/final_ann/cell_interaction/subset_total_500_rawt.h5","mat")
mat <- mydata$block0_values
rownames(mat) <- mydata$axis0
colnames(mat) <- mydata$axis1
mat <- Matrix(mat, sparse = TRUE)
meta <- read.csv("J:/scanpy_out/final_ann/cell_interaction/subset_total_500_meta.csv", row.names = 1)
sce <- CreateSeuratObject(mat,assay='RNA',meta.data=meta)
sce2 <- sce

# import snRNA-seq data
sce <- readRDS('J:/sn_data/sn_total_raw.RDS')
barcodes <- c()
for (i in  c('Neuron')) {
  print(i)
  sce1 <- subset(sce, subset = celltype==i)
  d <- 1:length(Cells(sce1))
  d1 <- sample(d, 700)
  sce1 <- sce1[,d1]
  barcodes <- c(barcodes, Cells(sce1))
}
sce1 <- sce[,barcodes]
sce1$ID <- 'Extra'
sce1$orig.ident <- sce1$nCount_RNA <- sce1$nFeature_RNA <- sce1$percent.mt <- sce1$celltype2 <- NULL
sce <- merge(sce1, sce2)
saveRDS(sce, "J:/scanpy_out/final_ann/cell_interaction/total_reference_raw2.RDS")

#--------------------------------------------------------------------------------------- Infering cell-cell interaction ----------------------------------------------------------------------------
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
#Convert("J:/scanpy_out/final_ann/cell_interaction/total_reference_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
sce.integrated <- readRDS("J:/scanpy_out/final_ann/cell_interaction/total_reference_raw2.RDS")
#Idents(sce.integrated) <- sce.integrated$celltype
sce.integrated$celltype <- factor_order_change(c('Telocyte_VSTM2A','Telocyte_NPY','FLC_ACTA2','FLC_SCARA5','FLC_CCL19','FLC_GREM1','FLC_APOD','FLC_CLDN1','FLC_KCNN3',
                                                 'Pericyte_RGS5','Pericyte_MYH11','Pericyte_CCL2','EC_Lym','EC_Venous','EC_Arterial',
                                                 'Epi_stem', 'TA','Enterocytes', 'Enterocytes_SLC26A3','Enterocytes_LYZ', 'Enterocytes_BEST4', 'Goblet_RETNLB', 'Goblet_ACHE',  'Tuft', 'Enteroendocrine',
                                                 'Glial','Neuron',
                                                 'GCB_MKI67','GCB','ProB','Brm_IGL','Brm_FCRL4','Brm_CD267','Brm','Plasmablast','Plasma',
                                                 'LTi','Th17','Th2','Treg_like','Tfh_like','CD4T_naive','CD8T_naive','CD8Tcm','CD8Tgd_ZNF683','CD8Tgd_CD39','CD8Trm','CD8Tem','CD8Tem_CX3CR1','T_cycling','gdT','NK',
                                                 'M_IL1B','M_C1Q','M_LYVE1','M_APOE','M_cycling','cDC2','cDC1','DC_LAMP3','pDC','Mast'), sce.integrated$celltype)
Idents(sce.integrated) <- sce.integrated$celltype
#DotPlot(sce.integrated, features = c('PTPRC','CD3D','CD19','IGHA1','CD14','MAP2','COL3A1','APOD','CLDN1','S100B','XBP1','RGS5','PECAM1','CD4'))
sce.integrated <- set_celltype(sce.integrated, new.cluster.ids = c('FLC','FLC','FLC','FLC',
                                                                   'FLC','FLC','FLC_APOD','FLC_CLDN1','FLC','Pericyte',
                                                                   'Pericyte','Pericyte','Endo','Endo','Endo','Epi',
                                                                   'Epi','Epi','Epi','Epi','Epi','Epi',
                                                                   'Epi','Epi','Epi','Glial','ENC',
                                                                   'B','B',
                                                                   'B','B','B','B','B','B',
                                                                   'B','ILC','T','T','T','T',
                                                                   'T','T','T','T','T','T',
                                                                   'T','T','T','T','ILC','Myeloid',
                                                                   'Myeloid','Myeloid','Myeloid','Myeloid','Myeloid','Myeloid',
                                                                   'Myeloid','Myeloid','Myeloid'))
sce.integrated$celltype <- Idents(sce.integrated)
sce.integrated <- subset(sce.integrated, subset = celltype%in%c('FLC_APOD','FLC_CLDN1','Glial','Endo','Pericyte','ENC'))
sce.integrated <- SCTransform(sce.integrated, verbose = FALSE)
sce.integrated <- RunPCA(sce.integrated, verbose = FALSE)
sce.integrated <- RunUMAP(sce.integrated, dims = 1:30, verbose = FALSE)
Idents(sce.integrated) <- sce.integrated$celltype
DimPlot(sce.integrated,label = T)

cellchat <- createCellChat(object = sce.integrated, group.by = "celltype", assay = "SCT")
#dplyr::glimpse(CellChatDB$interaction)
# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 6) # do parallel
# Preprocessing the expression data for cell-cell communication analysis
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat@idents <- factor(cellchat@idents)
cellchat <- projectData(cellchat, PPI.human)
#----------------------------------------------------- Part II: Inference of cell-cell communication network ------------------------------------
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
#cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat, 'J:/scanpy_out/final_ann/cell_interaction/GVB_neuron_cellchat7.RDS')
setwd("J:/scanpy_out/final_ann/cell_interaction/inter_GVB7/")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
levels(cellchat@idents)
vertex.receiver = c(1,2,3,4,5)
pathways.show.all <- cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  #print(i)
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 4, height = 8, units = 'in', dpi = 400)
}

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 30)
ht1 + ht2


