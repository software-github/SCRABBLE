# This is the main R file to generate the results related our manuscript
# figure 7 and the supplementary figures related to figure 7.
# Please contact Tao Peng: pengt@email.chop.edu if you have any questions 
# about the scripts or data

# load the libraries 
library(dplyr)
library(vsn)
library(Seurat)
library(scater)
library(edgeR)
library(gridExtra)
library(R.matlab)
library(cowplot)
library(biomaRt)
library(data.table)
library(lattice)
library(scImpute)
library(SCRABBLE)
library(VennDiagram)
library(Rtsne)
library(DT)
library(ggpubr)
library(ggsignif)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(refGenome)
library(pheatmap)
library(RColorBrewer)
library(dendsort)
library(entropy)
library(DrImpute)
library(splatter)
library(RColorBrewer)
library(mcriPalettes)
library(plotly)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(class)

source("analysis_library.R")

# ---------------------------------------------------------------------------------
# Fetal Brain
# ---------------------------------------------------------------------------------
# define the tissue names
tissue_name <- "FetalBrain"

# load the scRNAseq data
rnaseq_sc_mca <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"_rm.batch_dge.txt"))

# extract gene name informatioin
inter_gene <- rnaseq_sc_mca$V1

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# extract the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# extract the index of the cells of the specific tissue
index_cell <- cell_sc_mca$Tissue == "Fetal_Brain"

# extract the cell informatioin
cell_info <- cell_sc_mca$Cell.name[index_cell]

# extract the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# define the data.frame of the meta data
meta_data <- as.data.frame(meta_data)

# extract the cluster information
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]

# get the common cells 
common_cell <- intersect(cell_info,colnames(data_rnaseq))

# get the index of the common cells
index_select <- match(cell_info, colnames(data_rnaseq))

# get the data
data_select <- data_rnaseq[,..index_select]

# set up the row names
rownames(data_select) <- inter_gene

# set up the column name of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = paste0("data_sc_bulk/sc_",tissue_name,".RData"))

# ---------------------------------------------------------------------------------
# Small Intestine
# ---------------------------------------------------------------------------------

# define tissue name
tissue_name <- "SmallIntestine"

# load the scRNAseq data
rnaseq_sc_mca1 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"1_rm.batch_dge.txt"))

rnaseq_sc_mca2 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"2_rm.batch_dge.txt"))

rnaseq_sc_mca3 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"3_rm.batch_dge.txt"))

# get the common genes
inter_gene <- intersect(intersect(rnaseq_sc_mca1$V1,rnaseq_sc_mca2$V1),rnaseq_sc_mca3$V1)

# get the data of the 1st batch
index_1 <- match(inter_gene,rnaseq_sc_mca1$V1) 

rnaseq_sc_mca_1 <- rnaseq_sc_mca1[index_1,]

# get the data of the 2nd batch
index_2 <- match(inter_gene,rnaseq_sc_mca2$V1) 

rnaseq_sc_mca_2 <- rnaseq_sc_mca2[index_2,]

# get the data of the 3rd batch
index_3 <- match(inter_gene,rnaseq_sc_mca3$V1) 

rnaseq_sc_mca_3 <- rnaseq_sc_mca3[index_3,]

# combine the datasets
rnaseq_sc_mca <- cbind(rnaseq_sc_mca_1,rnaseq_sc_mca_2,rnaseq_sc_mca_3)

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# get the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# get the index of cells of the meta data
index_cell <- cell_sc_mca$Tissue == "Small-Intestine"

# get the cell information
cell_info <- cell_sc_mca$Cell.name[index_cell]

# get the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# get the meta data as a data.frame
meta_data <- as.data.frame(meta_data)

# get the cluster information
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]

# get the common cells 
common_cell <- intersect(cell_info,colnames(data_rnaseq))

# get the index of the common cells
index_select <- match(cell_info, colnames(data_rnaseq))

# get the data
data_select <- data_rnaseq[,..index_select]

# set up the row names of the data
rownames(data_select) <- inter_gene

# set up the column names of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = paste0("data_sc_bulk/sc_",tissue_name,".RData"))

# ---------------------------------------------------------------------------------
# Kidney
# ---------------------------------------------------------------------------------

# define tissue name
tissue_name <- "Kidney"

# load the scRNAseq data from two batches
rnaseq_sc_mca1 <- fread("data_sc_mca/rmbatch_dge/",tissue_name,"1_rm.batch_dge.txt")

rnaseq_sc_mca2 <- fread("data_sc_mca/rmbatch_dge/",tissue_name,"2_rm.batch_dge.txt")

# get the common genes
inter_gene <- intersect(rnaseq_sc_mca1$V1,rnaseq_sc_mca2$V1)

# get the data of the first batch
index_1 <- match(inter_gene,rnaseq_sc_mca1$V1) 

rnaseq_sc_mca_1 <- rnaseq_sc_mca1[index_1,]

# get the data of the second batch
index_2 <- match(inter_gene,rnaseq_sc_mca2$V1) 

rnaseq_sc_mca_2 <- rnaseq_sc_mca2[index_2,]

# combine all data
rnaseq_sc_mca <- cbind(rnaseq_sc_mca_1,rnaseq_sc_mca_2)

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# get the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# get the index of cells
index_cell <- cell_sc_mca$Tissue == tissue_name

cell_info <- cell_sc_mca$Cell.name[index_cell]

# get the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# get the data.frame of the meta data
meta_data <- as.data.frame(meta_data)

# get the cluster informaiton
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]


# get the cell index
index_select <- match(cell_info, colnames(data_rnaseq))

# extract the data
data_select <- as.data.frame(data_rnaseq[,..index_select])

# set up the row names of the data
rownames(data_select) <- inter_gene

# set up the row names of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = "data_sc_bulk/sc_",tissue_name,".RData")

# ---------------------------------------------------------------------------------
# Liver
# ---------------------------------------------------------------------------------

# define tissue name
tissue_name <- "Liver"

# load the scRNAseq data from two batches
rnaseq_sc_mca1 <- fread("data_sc_mca/rmbatch_dge/",tissue_name,"1_rm.batch_dge.txt")

rnaseq_sc_mca2 <- fread("data_sc_mca/rmbatch_dge/",tissue_name,"2_rm.batch_dge.txt")

# get the common genes
inter_gene <- intersect(rnaseq_sc_mca1$V1,rnaseq_sc_mca2$V1)

# get the data of the first batch
index_1 <- match(inter_gene,rnaseq_sc_mca1$V1) 

rnaseq_sc_mca_1 <- rnaseq_sc_mca1[index_1,]

# get the data of the second batch
index_2 <- match(inter_gene,rnaseq_sc_mca2$V1) 

rnaseq_sc_mca_2 <- rnaseq_sc_mca2[index_2,]

# combine all data
rnaseq_sc_mca <- cbind(rnaseq_sc_mca_1,rnaseq_sc_mca_2)

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# get the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# get the index of cells
index_cell <- cell_sc_mca$Tissue == tissue_name

cell_info <- cell_sc_mca$Cell.name[index_cell]

# get the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# get the data.frame of the meta data
meta_data <- as.data.frame(meta_data)

# get the cluster informaiton
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]


# get the cell index
index_select <- match(cell_info, colnames(data_rnaseq))

# extract the data
data_select <- as.data.frame(data_rnaseq[,..index_select])

# set up the row names of the data
rownames(data_select) <- inter_gene

# set up the row names of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = "data_sc_bulk/sc_",tissue_name,".RData")

# ---------------------------------------------------------------------------------
# Lung
# ---------------------------------------------------------------------------------

# define the tissue name
tissue_name <- "Lung"

# load all data
rnaseq_sc_mca1 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"1_rm.batch_dge.txt"))

rnaseq_sc_mca2 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"2_rm.batch_dge.txt"))

rnaseq_sc_mca3 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"3_rm.batch_dge.txt"))

# get the common genes
inter_gene <- intersect(intersect(rnaseq_sc_mca1$V1,rnaseq_sc_mca2$V1),rnaseq_sc_mca3$V1)

# get the data of the 1st batch
index_1 <- match(inter_gene,rnaseq_sc_mca1$V1) 

rnaseq_sc_mca_1 <- rnaseq_sc_mca1[index_1,]

# get the data of the 2nd batch
index_2 <- match(inter_gene,rnaseq_sc_mca2$V1) 

rnaseq_sc_mca_2 <- rnaseq_sc_mca2[index_2,]

# get the data of the 3rd batch
index_3 <- match(inter_gene,rnaseq_sc_mca3$V1) 

rnaseq_sc_mca_3 <- rnaseq_sc_mca3[index_3,]

# combine the datasets
rnaseq_sc_mca <- cbind(rnaseq_sc_mca_1,rnaseq_sc_mca_2,rnaseq_sc_mca_3)

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# get the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# get the index of cells of the meta data
index_cell <- cell_sc_mca$Tissue == tissue_name

# get the cell information
cell_info <- cell_sc_mca$Cell.name[index_cell]

# get the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# get the meta data as a data.frame
meta_data <- as.data.frame(meta_data)

# get the cluster information
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]

# get the common cells 
common_cell <- intersect(cell_info,colnames(data_rnaseq))

# get the index of the common cells
index_select <- match(cell_info, colnames(data_rnaseq))

# get the data
data_select <- data_rnaseq[,..index_select]

# set up the row names of the data
rownames(data_select) <- inter_gene

# set up the column names of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = paste0("data_sc_bulk/sc_",tissue_name,".RData"))

# ---------------------------------------------------------------------------------
# Placenta
# ---------------------------------------------------------------------------------

# load the scRNAseq data
rnaseq_sc_mca1 <- fread("data_sc_mca/rmbatch_dge/PlacentaE14.1_rm.batch_dge.txt")

rnaseq_sc_mca2 <- fread("data_sc_mca/rmbatch_dge/PlacentaE14.2_rm.batch_dge.txt")

# get the common genes
inter_gene <- intersect(rnaseq_sc_mca1$V1,rnaseq_sc_mca2$V1)

# get the data of the first batch
index_1 <- match(inter_gene,rnaseq_sc_mca1$V1) 

rnaseq_sc_mca_1 <- rnaseq_sc_mca1[index_1,]

# get the data of the second batch
index_2 <- match(inter_gene,rnaseq_sc_mca2$V1) 

rnaseq_sc_mca_2 <- rnaseq_sc_mca2[index_2,]

# combine all data
rnaseq_sc_mca <- cbind(rnaseq_sc_mca_1,rnaseq_sc_mca_2)

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# get the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# get the index of cells
index_cell <- cell_sc_mca$Tissue == "Placenta"

cell_info <- cell_sc_mca$Cell.name[index_cell]

# get the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# get the data.frame of the meta data
meta_data <- as.data.frame(meta_data)

# get the cluster informaiton
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]


# get the cell index
index_select <- match(cell_info, colnames(data_rnaseq))

# extract the data
data_select <- as.data.frame(data_rnaseq[,..index_select])

# set up the row names of the data
rownames(data_select) <- inter_gene

# set up the row names of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = "data_sc_bulk/sc_Placenta.RData")

# ---------------------------------------------------------------------------------
# Spleen
# ---------------------------------------------------------------------------------

# define the name of tissue
tissue_name <- "Spleen"

# load the data
rnaseq_sc_mca <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"_rm.batch_dge.txt"))

# extract gene name informatioin
inter_gene <- rnaseq_sc_mca$V1

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# extract the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# extract the index of the cells of the specific tissue
index_cell <- cell_sc_mca$Tissue == "Fetal_Brain"

# extract the cell informatioin
cell_info <- cell_sc_mca$Cell.name[index_cell]

# extract the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# define the data.frame of the meta data
meta_data <- as.data.frame(meta_data)

# extract the cluster information
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]

# get the common cells 
common_cell <- intersect(cell_info,colnames(data_rnaseq))

# get the index of the common cells
index_select <- match(cell_info, colnames(data_rnaseq))

# get the data
data_select <- data_rnaseq[,..index_select]

# set up the row names
rownames(data_select) <- inter_gene

# set up the column name of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = paste0("data_sc_bulk/sc_",tissue_name,".RData"))

# ---------------------------------------------------------------------------------
# FetalLiver
# ---------------------------------------------------------------------------------

# define the tissue name
tissue_name <- "FetalLiver"

# load the data
rnaseq_sc_mca1 <- fread(paste0("data_sc_mca/rmbatch_dge/",tissue_name,"E14.1_rm.batch_dge.txt"))

# extract gene name informatioin
inter_gene <- rnaseq_sc_mca$V1

# load cell information
cell_sc_mca <- fread("data_sc_mca/rmbatch_dge/MCA_CellAssignments.csv",header = TRUE)

# extract the data
data_rnaseq <- rnaseq_sc_mca[,-1]

# extract the index of the cells of the specific tissue
index_cell <- cell_sc_mca$Tissue == "Fetal_Brain"

# extract the cell informatioin
cell_info <- cell_sc_mca$Cell.name[index_cell]

# extract the meta data
meta_data <- cell_sc_mca[index_cell,-1]

# define the data.frame of the meta data
meta_data <- as.data.frame(meta_data)

# extract the cluster information
meta_data$Cluster <- t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]

# get the common cells 
common_cell <- intersect(cell_info,colnames(data_rnaseq))

# get the index of the common cells
index_select <- match(cell_info, colnames(data_rnaseq))

# get the data
data_select <- data_rnaseq[,..index_select]

# set up the row names
rownames(data_select) <- inter_gene

# set up the column name of the meta data
rownames(meta_data) <- colnames(data_select)

# save the data
save(data_select,meta_data,file = paste0("data_sc_bulk/sc_",tissue_name,".RData"))

# ---------------------------------------------------------------------------------
# Prepare the data used for the imputation
# ---------------------------------------------------------------------------------

# Define the tissue name and file name of the bulk RNAseq data
data_tissue_file <- read.table(file = "data/bulk_data_name.txt", sep = ",")

# The tissue name 
common_tissue <- c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
                   "Spleen", "Placenta", "FetalLiver", "Lung")

# get the data for the imputation
for(i in c(1:8)){
  
  get_data(common_tissue[i], data_tissue_file)
  
}

# ---------------------------------------------------------------------------------
# Run DrImpute
# ---------------------------------------------------------------------------------

# common tissue
common_tissue <- c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
                   "Spleen", "Placenta", "FetalLiver", "Lung")

# do the imputation and save the data in the folder "/imputation_data"
for(i in c(1:8)){
  
  tissue_name <- common_tissue[i]
  
  data <- readRDS(paste0("/data_sc_bulk/",tissue_name,"_imputation.rds"))
  
  extdata <- DrImpute(data[[1]]) 
  
  saveRDS(extdata, file = paste0("/imputation_data/",
                                 tissue_name,
                                 "_drimpute_imputation.rds"))
  
}


# ---------------------------------------------------------------------------------
# Run scImpute
# ---------------------------------------------------------------------------------

# common tissue
common_tissue <- c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
                   "Spleen", "Placenta", "FetalLiver", "Lung")

# do the imputation and save the data in the folder "/imputation_data"
for(i in c(1:8)){
  
  tissue_name <- common_tissue[i]
  
  data <- readRDS(paste0("/data_sc_bulk/",tissue_name,"_imputation.rds"))
  
  write.table(as.matrix(data[[1]]), paste0("/imputation_data/",tissue_name,"_scimpute.csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  scimpute(
    paste0("/imputation_data/",tissue_name,"_scimpute.csv"),
    infile = "csv",           
    outfile = "csv",          
    out_dir = paste0("/imputation_data/scImpute_", tissue_name, "_"),
    drop_thre = 0.5,          
    Kcluster = 2,
    ncores = 2)             
  
  data_scimpute <- read.table( file = paste0("/imputation_data/scImpute_", tissue_name,
                                            "_scimpute_count.csv") ,
                              header = TRUE, sep=",")
  
  data_scimpute$X <- NULL
  
  system(paste0("rm -R imputation_data/scImpute_", tissue_name,"*"))
  
  # save the data
  saveRDS(data_scimpute, file = paste0("/imputation_data/",
                                      tissue_name,
                                      "_scimpute_imputation.rds"))
  
}


# ---------------------------------------------------------------------------------
# Run MAGIC
# ---------------------------------------------------------------------------------

# import magic
# import pandas as pd
# import matplotlib.pyplot as plt
# import os
# 
# cwd = os.getwd()
# 
# def run_magic_singlecell(tissue_name):
#   filename = '/imputation_data/'+tissue_name+'_scimpute.csv'
#   X = pd.read_csv(filename)
#   magic_operator = magic.MAGIC()
#   X_magic = magic_operator.fit_transform(X.T)
#   out_magic = X_magic.T
#   out_magic.to_csv(cwd+'/imputation_data/'+tissue_name+'_magic_imputation.csv', sep = '\t', header= None)
# 
# 
# for i in range(0,8):
#   X = pd.read_csv('common_tissue.csv')
#   tissue_name_all = X.x
#   run_magic_singlecell(tissue_name_all[i]) 


# ---------------------------------------------------------------------------------
# Run SCRABBLE
# ---------------------------------------------------------------------------------

# common tissue
common_tissue <- c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
                   "Spleen", "Placenta", "FetalLiver", "Lung")

# parameter settings
parameterT <- rbind(c(1,1e-7,1e-1),
                    c(1,1e-7,1e0),
                    c(1,1e-7,1e0),
                    c(1,1e-7,1e-1),
                    c(1,1e-7,1e0),
                    c(1,1e-7,1e-1),
                    c(1,1e-7,1e-1),
                    c(1,1e-7,1e-2))

# do the imputation and save the data in the folder "/imputation_data
for(i in c(1:8)){
  
  tissue_select <- common_tissue[i]
  
  data <- readRDS(paste0("/data_sc_bulk/",tissue_select,"_imputation.rds"))
  
  parameter <- parameterT[i,]
  
  result <- scrabble(data,
                     parameter = parameter, 
                     nIter = 30,
                     error_out_threshold = 1e-7, 
                     nIter_inner = 30,
                     error_inner_threshold = 1e-5)
  
  saveRDS(result, 
          file = paste0("/imputation_data/",
                        tissue_select,
                        "_scrabble_imputation_new.rds"))
  
}


# ---------------------------------------------------------------------------------
# Run tSNE
# ---------------------------------------------------------------------------------

# common tissue
common_tissue <- c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
                   "Spleen", "Placenta", "FetalLiver", "Lung")

# tSNE parameters
initial_dims_value <- 60
perplexity_value <- 100

for(i in c(1:8)){
  
  tissue_name <- common_tissue[i]
  
  data_tsne <- list()
  
  # load the raw data and calculate the tSNE matrix 
  data <- readRDS(paste0(cwd,"/data_sc_bulk/",
                         tissue_name,
                         "_imputation.rds"))
  
  set.seed(42)
  
  data_tsne[[1]] <- Rtsne(t(as.matrix(data[[1]])), 
                          initial_dims = initial_dims_value, 
                          perplexity = perplexity_value)
  
  
  # load the DrImptue data and calculate the tSNE matrix
  data_drimpute <- readRDS(file = paste0("imputation_data/",
                                         tissue_name,
                                         "_drimpute_imputation.rds"))
  
  set.seed(42)
  
  data_tsne[[2]] <- Rtsne(t(as.matrix(data_drimpute)), 
                          initial_dims = initial_dims_value, 
                          perplexity = perplexity_value)
  
  # load the scImpute data and calculate the tSNE matrix
  data_scimpute <- readRDS(file = paste0("imputation_data/",
                                         tissue_name,
                                         "_scimpute_imputation.rds"))
  
  set.seed(42)
  
  data_tsne[[3]] <- Rtsne(t(as.matrix(data_scimpute)), 
                          initial_dims = initial_dims_value, 
                          perplexity = perplexity_value)
  
  
  # load the imputation results by MAGIC and calculate the tSNE matrix
  data_magic <- read.table( file = paste0("imputation_data/",
                                          tissue_name,
                                          "_magic_imputation.csv"),
                            header = FALSE, sep="\t")
  
  data_magic <- data_magic[,-1]
  
  set.seed(42)
  
  data_tsne[[4]] <- Rtsne(t(as.matrix(data_magic)), 
                          initial_dims = initial_dims_value, 
                          perplexity = perplexity_value)
  
  # load the SCRABBLE data and calculate the tSNE matrix
  data_scrabble <- readRDS(file = paste0("imputation_data/",
                                         tissue_name,
                                         "_scrabble_imputation.rds"))
  
  set.seed(42)
  
  data_tsne[[5]] <- Rtsne(t(as.matrix(data_scrabble)), 
                          initial_dims = initial_dims_value, 
                          perplexity = perplexity_value)
  
  saveRDS(data_tsne, file = paste0("data_all/",tissue_name,"_TSNE.rds"))
  
}

# ---------------------------------------------------------------------------------
# Plot tSNE and Bun Index
# ---------------------------------------------------------------------------------

tissue_all <- c("FetalBrain",  "FetalLiver")

for(i in c(1:2)){

  pdf_dun_tsne1(tissue_all[i])
  
}


tissue_all <- c("Kidney",  "Lung", "Placenta", "SmallIntestine", "Liver", "Spleen")
for(i in c(1:6)){
  
  pdf_dun_tsne(tissue_all[i])
  
}

# ---------------------------------------------------------------------------------










