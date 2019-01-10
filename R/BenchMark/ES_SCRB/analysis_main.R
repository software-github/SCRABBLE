# This is the main R file to generate the results related our manuscript
# figure 4AB and the supplementary figures related to figure 4AB.
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

# Load the data
# ---------------------------------------------------------------------------------
# load the single cell RNAseq data, bulk RNAseq data, and the gene name information
data_sc_info <- read.table(file = "data/data_sc.txt", 
                           fill = TRUE, 
                           header = TRUE, 
                           sep = "\t", 
                           stringsAsFactors = FALSE)

# load the annotationn file
anno_temp <- read.table(file = "data/mouse_gene_id.txt", 
                        header = FALSE, 
                        sep = ",", 
                        stringsAsFactors = FALSE)

anno <- data.frame(apply(anno_temp,2,function(x)gsub('\\s+', '',x)),stringsAsFactors = FALSE)

# load the bulk data
data_bulk_info <- read.table(file = "data/data_bulk.txt", 
                             fill = TRUE, 
                             header = TRUE, 
                             sep = "\t", 
                             stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------

# Preprocess the data
# ---------------------------------------------------------------------------------
# prepaare the data
data_sc_info_tmp <- cbind(id = rownames(data_sc_info),data_sc_info)

# make the annotation of genes
index <- match(rownames(data_sc_info),anno$V1)

gene_name <- anno[index,]

index <- match(gene_name$V1,rownames(data_sc_info))

# make the gene names
data_sc_tmp <- data_sc_info_tmp[index,]

data_sc_tmp <- cbind(gene_name,data_sc_tmp)

data_sc <- data_sc_tmp[rowSums(is.na(data_sc_tmp)) == 0,]

# extract the common genes
index <- match(data_bulk_infos$tracking_id, data_sc$V2)

gene_common <- data_sc$V2[index]

gene_common <- gene_common[!is.na(gene_common)]

# select the scRNAseq data
index <- match(gene_common,data_sc$V2)

data_sc_final <- data_sc[index,]

# select the bulk RNAseq data
index <- match(gene_common,data_info1$tracking_id)

data_bulk_final <- data_info1[index,]

# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# select the bulk RNAseq, SCRB RNAseq data, Dropseq data, and gene names

# select the bulk RNAseq data
data_bulk_rep <- data.frame(rep1 = data_bulk_final$WT_ES, stringsAsFactors = FALSE)

# select the DropSeq scRNAseq data
col_name <- colnames(data_sc_final)
is.Dropseq <- grepl("^Drop", col_name)
data_sc_dropseq <- data_sc_final[,is.Dropseq]

# select the SCRB scRNAseq data
is.SCRB <- grepl("^SCRB", col_name)
data_sc_SCRBseq <- data_sc_final[,is.SCRB]

# select the gene names
gene_bulk_sc <- data_sc_final$V2

# ---------------------------------------------------------------------------------

# normalize by the library size
# ---------------------------------------------------------------------------------
# by the bulk RNAseq data
tmp <-  rep.row(colSums(data_bulk_rep),dim(data_bulk_rep)[1])
data_bulk_rep_normalized <- 1e6*data_bulk_rep/tmp

# by the dropseq RNAseq data
tmp <-  rep.row(colSums(data_sc_dropseq),dim(data_sc_dropseq)[1])
data_sc_dropseq_normalized <- 1e6*data_sc_dropseq/tmp

# by the scrb RNAseq data
tmp <-  rep.row(colSums(data_sc_SCRBseq),dim(data_sc_SCRBseq)[1])
data_sc_SCRBseq_normalized <- 1e6*data_sc_SCRBseq/tmp

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# calculate the number cells of Dropseq RNAseq data and of the SCRB RNAseq data
n_drop <- dim(data_sc_dropseq_normalized)[2]
n_SCRB <- dim(data_sc_SCRBseq_normalized)[2]

# get the index of genes of Dropseq RNAseq data and of the SCRB RNAseq data without dropouts
tmp_drop <- rowSums((data_sc_dropseq_normalized) > 0) == n_drop
tmp_SCRB <- rowSums((data_sc_SCRBseq_normalized) > 0) == n_SCRB

# calculate the mean values of bulk RNAseq, Dropseq RNAseq, and SCRB RNAseq
means_drop <- rowMeans(data_sc_dropseq_normalized)
means_SCRB <- rowMeans(data_sc_SCRBseq_normalized)
means_bulk <- rowMeans(data_bulk_rep_normalized)

# get the common gene between the dropseq and scrb
tmp_index <- tmp_drop & tmp_SCRB
common_gene_single <- gene_bulk_sc[tmp_index]

# extract the common genes without dropouts only at SCRB RNAseq data
gene_only_SCRB <- gene_bulk_sc[tmp_SCRB]



# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# calculate the normalized ratio between Dropseq and SCRB
x <- means_drop[tmp_index]
y <- means_SCRB[tmp_index]
c1 <- lm(x ~ y + 0)

# calculate the normalized ratio between Dropseq and bulk
x <- means_drop[tmp_index]
y <- means_bulk[tmp_index]
c2 <- lm(x ~ y + 0)

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# normalize the data
data_sc_SCRBseq_normalized <- c1$coefficients*data_sc_SCRBseq_normalized 
data_bulk_rep_normalized <- c2$coefficients*data_bulk_rep_normalized

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# remove the genes with low expression values
rowMeans1 <- rowMeans(data_bulk_rep_normalized)

index <- rowMeans1 > 0.5

data_sc_dropseq_tmp <- data_sc_dropseq_normalized[index,]

data_sc_SCRB_tmp <- data_sc_SCRBseq_normalized[index,]

gene_list <- gene_bulk_sc[index]

data_bulk_tmp <- rowMeans1[index]

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# remove the mt genes
is.mito <- grepl("^mt", gene_list)

gene_list <- gene_list[which(!is.mito)]

data_sc_dropseq_tmp <- data_sc_dropseq_tmp[which(!is.mito),]

data_sc_SCRB_tmp <- data_sc_SCRB_tmp[which(!is.mito),]

data_bulk_tmp <- data_bulk_tmp[which(!is.mito)]

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# calculate the variance data of dropseq data
tmp <- rowVars(as.matrix(data_sc_dropseq_tmp), suma = NULL, std = FALSE)

tmp[is.na(tmp)] <- 0

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# remove the genes with low variances
index <- which(tmp > 1e-2)

gene_list <- gene_list[index]

data_sc_dropseq <- data_sc_dropseq_tmp[index,]

data_sc_scrb <- data_sc_scrb_tmp[index,]

data_bulk <- data_bulk_tmp[index]

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# save the data
saveRDS(gene_list, file = "/data_all/impute_gene.rds")

saveRDS(gene_only_SCRB, file = "/data_all/gene_only_SCRB.rds")

write.table(data_sc_dropseq, 
            file = paste0("/data_all/sc_data_ES.csv"),
            sep=',',row.names = FALSE, col.names = FALSE)

write.table(data_sc_scrb, 
            file = paste0("/data_all/sc_data_ES_SCRB.csv"),
            sep=',',row.names = FALSE, col.names = FALSE)

write.table(data_bulk, 
            file = paste0("/data_all/bulk_data_ES.csv"),
            sep=',',row.names = FALSE, col.names = FALSE)

# ---------------------------------------------------------------------------------
# Run DrImpute 
# ---------------------------------------------------------------------------------

dir.create(file.path("imputation_data/"), 
           showWarnings = FALSE)

data_sc <- as.matrix(fread("/data_all/sc_data_ES.csv"))

extdata <- DrImpute(data_sc)

saveRDS(extdata, file = "imputation_drimpute_data/data_drimpute_imputation.rds")

# ---------------------------------------------------------------------------------
# Run scImpute 
# ---------------------------------------------------------------------------------

# create the folder of the imputation result
dir.create(file.path("imputation_data/"), 
           showWarnings = FALSE)

# load the data
data_sc <- as.matrix(fread("/data_all/sc_data_ES.csv"))

write.table(data_sc, file = "/data_all/sc_data_ES_scimpute.csv", sep = ",")

scimpute(
  "/data_all/sc_data_ES_scimpute.csv", 
  infile = "csv",           
  outfile = "csv",         
  out_dir = paste0("imputation_data/"),          
  drop_thre = 0.5,          
  Kcluster = 2,
  ncores = 2)             

# move all data to the  data folder
tmp1 <- read.table( file = "scimpute_count.csv" , header = TRUE, sep=",")
tmp1$X <- NULL
tmp1[is.na(tmp1)] <- 0

# save the data
write.table(tmp1, file = "imputation_data/data_scimpute_imputation.csv",
            sep=',',row.names = FALSE, col.names = FALSE)


# ---------------------------------------------------------------------------------
# Run MAGIC
# ---------------------------------------------------------------------------------

# cwd = os.getwd()
# 
# if not os.path.exists(cwd+"/imputation_data"):
#   os.makedirs(cwd+"/imputation_data")
# 
# 
# X = pd.read_csv(cwd + "/data_all/sc_data_ES.csv",sep = ',',header=None)
# 
# magic_operator = magic.MAGIC(n_pca = 76)
# 
# X_magic = magic_operator.fit_transform(X.T)
# 
# out_magic = X_magic.T
# 
# out_magic.to_csv(cwd+"/magic_data/data_imputation.csv", sep = ',', header= None)


# ---------------------------------------------------------------------------------
# Run SCRABBLE
# ---------------------------------------------------------------------------------

# create the folder
dir.create(file.path("/imputation_data/"), 
           showWarnings = FALSE)

# load the scRNAseq data
data_sc <- as.matrix(fread(paste0(cwd,"/data_all/sc_data_ES.csv")))

# load the bulk RNAseq data
data_bulk <- as.matrix(fread(paste0(cwd,"/data_all/bulk_data_ES.csv")))

# prepare the data 
data <- list()
data[[1]] <- as.matrix(data_sc)
data[[2]] <- as.matrix(data_bulk)

# set up the parameter
parameter <- c(1,1e-7,1e-5)

# run the SCRABBLE
result <- scrabble(data,
                   parameter = parameter, 
                   nIter = 30,
                   error_out_threshold = 1e-14, 
                   nIter_inner = 30,
                   error_inner_threshold = 1e-14)

# save the result
saveRDS(result, file = "/imputation_data/data_scrabble_imputation.rds")

# ---------------------------------------------------------------------------------
# Assemble all data: raw data, gold standard data, imputed data
# ---------------------------------------------------------------------------------
# load DropSeq data
data_sc_Dropseq <- read.table(file = "data_all/sc_data_ES.csv", 
                              fill = TRUE, header = FALSE, sep = ",", 
                              stringsAsFactors = FALSE)
data_sc_Dropseq <- log10(as.matrix(data_sc_Dropseq) + 1)

# load SCRB RNAseq data
data_sc_SCRBseq <- read.table(file = "data_all/sc_data_ES_SCRB.csv", 
                              fill = TRUE, header = FALSE, sep = ",", 
                              stringsAsFactors = FALSE)
data_sc_SCRBseq <- log10(as.matrix(data_sc_SCRBseq) + 1)

# load bulk RNAseq data
data_sc_bulk <- read.table(file = "data_all/bulk_data_ES.csv", 
                           fill = TRUE, header = FALSE, sep = ",", 
                           stringsAsFactors = FALSE)
data_sc_bulk <- log10(as.matrix(data_sc_bulk) + 1)

# load the imputed data of DrImpute
data_sc_drimpute <- readRDS(file = "imputation_data/data_drimpute_imputation.rds")

data_sc_drimpute <- log10(data_sc_drimpute + 1)

# load the imputed data of scImpute
data_sc_scimpute <- read.table(file = "imputation_data/data_scimpute_imputation.csv", 
                               fill = TRUE, header = FALSE, sep = ",", 
                               stringsAsFactors = FALSE)

data_sc_scimpute <- log10(data_sc_scimpute + 1)

# load the imputed data of MAGIC
data_sc_magic <- read.table(file = "imputation_data/data_magic_imputation.csv", 
                            header = FALSE, sep = ",", 
                            stringsAsFactors = FALSE)

data_sc_magic <- data_sc_magic[,-1]

data_sc_magic <- log10(data_sc_magic + 1)

# load the imputed data of SCRABBLE
data_sc_scrabble <- readRDS(file = "imputation_data/data_scrabble_imputation.rds")

data_sc_scrabble <- log10(data_sc_scrabble + 1) 

# load the gene lists
gene_filter <- readRDS(file = "/data_all/impute_gene.rds")

gene_only_SCRB <- readRDS(file = "/data_all/gene_only_SCRB.rds")

# combine all data
data_de <- cbind(data_sc_SCRBseq,
                 data_sc_Dropseq,
                 data_sc_drimpute,
                 data_sc_scimpute,
                 data_sc_magic,
                 data_sc_scrabble)

# define cell types
dataType <- c(rep(1,dim(data_sc_Dropseq)[2]),
              rep(2,dim(data_sc_SCRBseq)[2]),
              rep(3,dim(data_sc_Dropseq)[2]),
              rep(4,dim(data_sc_Dropseq)[2]),
              rep(5,dim(data_sc_Dropseq)[2]),
              rep(6,dim(data_sc_Dropseq)[2]))

# ---------------------------------------------------------------------------------
# Plot the boxplot for the comparison
# ---------------------------------------------------------------------------------

p <- performance_comparison(gene_only_SCRB, gene_filter, 22, data_de)

ggsave(
  filename = paste0("Figure_performance_59_dropout_29.pdf"),
  plot = p,
  width = 4,
  height = 4
)

p <- performance_comparison(gene_only_SCRB, gene_filter, 30, data_de)

ggsave(
  filename = paste0("Figure_performance_59_dropout_29.pdf"),
  plot = p,
  width = 4,
  height = 4
)

# ---------------------------------------------------------------------------------
# Plot the distribution comparison
# ---------------------------------------------------------------------------------
# plot the distribution
p <- performance_distribution_comparison(gene_only_SCRB, gene_filter, 22, data_de)

# Save the PDF of the plot of the distribution of the first 20 genes
main <- grid.arrange(grobs = p[c(1:20)],ncol = 4, top = "main")

ggsave(
  filename = paste0("Figure_gene_performance_compare_1.pdf"),
  plot = main,
  width = 14,
  height = 20
)

# Save the PDF of the plot of the distribution of the next 20 genes
main <- grid.arrange(grobs = p[c(21:40)],ncol = 4, top = "main")

ggsave(
  filename = paste0("Figure_gene_performance_compare_2.pdf"),
  plot = main,
  width = 14,
  height = 20
)

# Save the PDF of the plot of the distribution of the final 16 genes
main <- grid.arrange(grobs = p[c(41:56)],ncol = 4, top = "main")

ggsave(
  filename = paste0("Figure_gene_performance_compare_3.pdf"),
  plot = main,
  width = 14,
  height = 20
)












