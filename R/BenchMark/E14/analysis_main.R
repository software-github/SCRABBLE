# This is the main R file to generate the results related our manuscript
# figure 4CD and the supplementary figures related to figure 4CD.
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
# load the bulk RNAseq data
data_bulk0 <- read.table(file = "data/bulkRNAseq.csv", 
                         header = F, sep = ",",stringsAsFactors = F)

# load the scRNAseq data
data_sc0 <- data.frame(read.table(file = "data/scRNAseq.csv", 
                                  header = F, sep = ",",stringsAsFactors = F))

# get the gene names of bulk RNAseq and scRNAseq data
common_index <- which(data_bulk0$V1 != data_sc0$V1)

common_gene <- data.frame(gene = data_bulk0$V1[common_index],stringsAsFactors = F)

df <- transform(common_gene, test=do.call(rbind, strsplit(gene, '-', fixed=TRUE)), stringsAsFactors=F)

df$gene1 <- paste0(as.numeric(df$test.2),'-',df$test.1)

data_sc0$V1[common_index] <- df$gene

data_bulk1 <- data.frame(rep1 = data_bulk0$V2, rep2 = data_bulk0$V3, stringsAsFactors = F)

data_sc1 <- data_sc0[,-1]

gene <- make.unique(data_bulk0$V1)

gene[which(gene == 'Zfp42')] <- "Rex1"
gene[which(gene == 'Prdm1')] <- "Blimp1"

# remove the low quality of cells
sumV <- colSums(data_sc1)

data_sc2 <- data_sc1[,sumV > 2000]

# normalize the data by the library size
tmp <-  rep.row(colSums(data_bulk1),dim(data_bulk1)[1])

data_bulk_normalized <- 1e6*data_bulk1/tmp

tmp <-  rep.row(colSums(data_sc2),dim(data_sc2)[1])

data_sc_normalized <- 1e6*data_sc2/tmp

# select the genes without dropouts
n_sc <- dim(data_sc_normalized)[2]

tmp_sc <- rowSums((data_sc_normalized) > 0) == n_sc

means_sc <- rowMeans(data_sc_normalized)

means_bulk <- rowMeans(data_bulk_normalized)

# calculate the ratio between scRNAseq and bulk RNAseq
x <- means_sc[tmp_sc]

y <- means_bulk[tmp_sc]

c1 <- lm(x ~ y + 0)

# calculate the normalized data
data_sc3 <- data_sc_normalized

data_bulk2 <- c1$coefficients*data_bulk_normalized

# save all the data
write.table(data_sc3,"data_all/data_sc.csv",
            sep=',',row.names = F, col.names = F )

write.table(data_bulk2,"data_all/data_bulk.csv",
            sep=',',row.names = F, col.names = F )

write.table(gene,"data_all/gene.csv",
            sep=',',row.names = F, col.names = F )

# ---------------------------------------------------------------------------------
# Run DrImpute 
# ---------------------------------------------------------------------------------

dir.create(file.path("imputation_drimpute_data/"))

data_sc <- fread("/data_all/data_sc.csv")

extdata <- DrImpute(log10(as.matrix(data_sc) + 1))

extdata <- 10^extdata - 1

saveRDS(extdata, file = "imputation_drimpute_data/data_drimpute_imputation.rds")

# ---------------------------------------------------------------------------------
# Run scImpute 
# ---------------------------------------------------------------------------------

# create the folder of the imputation result
dir.create(file.path("imputation_scimpute_data/"))

# load the data
data_sc <- as.matrix(fread("/data_all/data_sc.csv"))

write.table(data_sc, file = "/data_all/sc_data_scimpute.csv", sep = ",")

scimpute(
  "/data_all/sc_data_scimpute.csv", 
  infile = "csv",           
  outfile = "csv",         
  out_dir = paste0(cwd1),          
  drop_thre = 0.5,          
  Kcluster = 2,
  ncores = 2)             

# move all data to the  data folder
tmp1 <- read.table( file = "scimpute_count.csv" , header = TRUE, sep=",")
tmp1$X <- NULL
tmp1[is.na(tmp1)] <- 0

# save the data
write.table(tmp1, file = "imputation_drimpute_data/data_scimpute_imputation.csv",
            sep=',',row.names = FALSE, col.names = FALSE)



# ---------------------------------------------------------------------------------
# Run MAGIC
# ---------------------------------------------------------------------------------

# cwd = os.getwd()
# 
# if not os.path.exists(cwd+"/magic_data"):
#   os.makedirs(cwd+"/magic_data")
# 
# 
# X = pd.read_csv(cwd + "/data_all/data_sc.csv",sep = ',',header=None)
# 
# magic_operator = magic.MAGIC()
# 
# X_magic = magic_operator.fit_transform(X.T)
# 
# out_magic = X_magic.T
# 
# out_magic.to_csv(cwd+"/magic_data/data_magic_imputation.csv", sep = ',', header= None)


# ---------------------------------------------------------------------------------
# Run SCRABBLE
# ---------------------------------------------------------------------------------

# create the folder
dir.create(file.path("/imputation_scrabble_data/"))

# load the scRNAseq data
data_sc <- as.matrix(fread(paste0(cwd,"/data_all/data_sc.csv")))

# load the bulk RNAseq data
data_bulk <- as.matrix(fread(paste0(cwd,"/data_all/data_bulk.csv")))

# prepare the data 
data <- list()
data[[1]] <- as.matrix(data_sc)
data[[2]] <- as.matrix(data_bulk)

# set up the parameter
parameter <- c(1,1e-7,1e-3)

# run the SCRABBLE
result <- scrabble(data,
                   parameter = parameter, 
                   nIter = 30,
                   error_out_threshold = 1e-14, 
                   nIter_inner = 30,
                   error_inner_threshold = 1e-14)

# save the result
saveRDS(result, file = "/imputation_scrabble_data/data_scrabble_imputation.rds")

# ---------------------------------------------------------------------------------
# Assemble all data: raw data, gold standard data, imputed data
# ---------------------------------------------------------------------------------

data_sc <- read.table( file = paste0("data_all/data_sc.csv") , 
                       header = F, sep=",", stringsAsFactors = F)

data_bulk <- read.table(file = paste0("data_all/data_bulk.csv"), 
                        fill = TRUE, header = FALSE, sep = ",", stringsAsFactors = FALSE)

data_sc_drimpute <- readRDS(file = "imputation_drimpute_data/data_drimpute_imputation.rds")

data_scrabble <- readRDS(file = "imputation_scrabble_data/data_scrabble_imputation.rds")


data_sc_magic <- read.table(file = paste0("imputation_magic_data/data_magic_imputation.csv"), 
                            header = FALSE, sep = ",", stringsAsFactors = FALSE)

data_sc_magic <- data_sc_magic[,-1]

data_sc_magic[data_sc_magic < 0] <- 0


data_sc_scimpute <- read.table(file = "imputation_scimpute_data/data_scimpute_imputation.csv", 
                               fill = TRUE, header = FALSE, sep = ",", stringsAsFactors = FALSE)

gene_name <- read.table( file = "data_all/gene_impute.csv", 
                         header = F, sep=",", stringsAsFactors = F)

gene_FISH <- c("Blimp1","Esrrb","Tbp","Tbx3","Klf4","Pecam1","Prdm14","Sox2","Nanog","Tet1","Nr0b1","Tbp")


# ---------------------------------------------------------------------------------
# Normalize the smFISH data
# ---------------------------------------------------------------------------------
# Sox2 is used to normaliz the data
gene1 <- "Sox2"

tmp <- data.frame(t(data_sc[which(gene_name$V1 == gene1),]))

colnames(tmp) <- "value"

tmp <- tmp[tmp$value > 0,]

tmp <- log10(tmp$value + 1)

data_FISH <- read.table( file = paste0("data_mFISH/",gene1,".csv") , 
                         header = F, sep=",", stringsAsFactors = F)

data_FISH <- log10(data_FISH + 1)

ratio <- mean(tmp)/mean(data_FISH$V1)

# ---------------------------------------------------------------------------------
# Normalize the smFISH data
# ---------------------------------------------------------------------------------
# define the KS statistics
ks1 <- c()
ks2 <- c()
ks3 <- c()
ks4 <- c()
ks5 <- c()
ks6 <- c()

for (i in c(1:length(gene_FISH))){
  
  # get the index of gene
  gene1 <- gene_FISH[i]
  
  index <- which(gene_name$V1 == gene1)
  
  # get the raw data
  tmp <- data.frame(t(data_sc[index,]))
  
  colnames(tmp) <- "value"
  
  tmp1 <- log10(tmp+1)
  
  tmp1$e <- 'Raw'
  
  # get the imputed data of DrImpute
  tmp <- data.frame(data_sc_drimpute[index,])
  
  colnames(tmp) <- "value"
  
  tmp2 <- log10(tmp+1)
  
  tmp2$e <- "DrImpute"
  
  # get the imputed data of scImpute
  tmp <- data.frame(t(data_sc_scimpute[index,]))
  
  colnames(tmp) <- "value"
  
  tmp3 <- log10(tmp+1)
  
  tmp3$e <- "scImpute"
  
  # get the imputed data of magic
  tmp <- data.frame(t(data_sc_magic[index,]))
  
  colnames(tmp) <- "value"
  
  tmp4 <- log10(tmp+1)
  
  tmp4$e <- "MAGIC"
  
  # get the imputed data of scrabble
  tmp <- data.frame((data_sc_scrabble[index,]))
  
  colnames(tmp) <- "value"
  
  tmp5 <- log10(tmp+1)
  
  tmp5$e <- "SCRABBLE"
  
  # get the FISH data
  data_FISH <- read.table( file = paste0("data_mFISH/",gene1,".csv"), 
                           header = F, sep=",", stringsAsFactors = F)
  
  tmp <- ratio*log10(data_FISH+1)
  
  colnames(tmp) <- "value"
  
  tmp6 <- tmp
  
  tmp6$e <- "smFISH"
  
  # calculate the KS statistics
  # Raw vs smFISH
  a <- ks.test(tmp1$value, tmp6$value)
  
  ks1 <- rbind(ks1,a[[1]])
  
  # drimpute vs smFISh
  a <- ks.test(tmp2$value, tmp6$value)
  
  ks2 <- rbind(ks2,a[[1]])
  
  # scimpute vs smFISH
  a <- ks.test(tmp3$value, tmp6$value)
  
  ks3 <- rbind(ks3,a[[1]])
  
  # magic vs smFISH
  a <- ks.test(tmp4$value, tmp6$value)
  
  ks4 <- rbind(ks4,a[[1]])
  
  # scrabble vs smFISH
  a <- ks.test(tmp5$value, tmp6$value)
  
  ks5 <- rbind(ks5,a[[1]])
  
}

# assemble KS statistics between raw and smFISH
ks11 <- data.frame(ks = unlist(ks1))

ks11$group <- 1

# assemble KS statistics between drimpute and smFISH
ks22 <- data.frame(ks = unlist(ks2))

ks22$group <- 2

# assemble KS statistics between scimpute and smFISH
ks33 <- data.frame(ks = unlist(ks3))

ks33$group <- 3

# assemble KS statistics between magic and smFISH
ks44 <- data.frame(ks = unlist(ks4))

ks44$group <- 4

# assemble KS statistics between scrabble and smFISH
ks55 <- data.frame(ks = unlist(ks5))

ks55$group <- 5

# combine all KS statistics
ksData <- data.frame(rbind(as.matrix(ks11),
                           as.matrix(ks22),
                           as.matrix(ks33),
                           as.matrix(ks44),
                           as.matrix(ks55)))

colnames(ksData) <- c("ks","group")

# set up the comparison
my_comparisons <- list( c("1", "5"),c("2", "5"), c("3", "5"), c("4", "5"))

# calculate pvalues
pval1 <- compare_means(ks ~ group,data = ksData, method = "t.test", paired = T)

pval <- c(paste0('p1 = ',formatC(pval1$p[4], format = "e", digits = 2)),
          paste0('p1 = ',formatC(pval1$p[7], format = "e", digits = 2)),
          paste0('p1 = ',formatC(pval1$p[9], format = "e", digits = 2)),
          paste0('p1 = ',formatC(pval1$p[10], format = "e", digits = 2)))

# plot the boxplot
pp <- ggboxplot(ksData, x = "group", y = "ks", fill = "group",
                palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00")) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  ylim(c(0,1.5)) + 
  theme_bw() +
  geom_signif(comparisons = my_comparisons, 
              annotations = pval,
              tip_length = 0.03,
              y_position = c(1.4,1.3, 1.2, 1.1)) +
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  xlab("Data Sets") + 
  ylab("K-S Statistics")

# save the PDF file of the boxplot
ggsave(filename=paste0("Figure_smFISH_boxplot.pdf"), 
       plot=pp, 
       width = 4, 
       height = 3)

# ---------------------------------------------------------------------------------
# Plot the distributions of comparison
# ---------------------------------------------------------------------------------
# define the plot list
p <- list()

for (i in c(1:length(gene2))){
  
  gene1 <- gene2[i]
  
  p[[i]] <- distribution_FISH(data_sc, 
                              data_drimpute, 
                              data_scimpute, 
                              data_magic, 
                              data_scrabble, 
                              gene_name, 
                              gene1,
                              ratio) 
}

# assemble the graph
main <- grid.arrange(grobs = p,ncol = 4, top = "main")

# save PDF file of distribution plots
ggsave(filename=paste0("Figure_smFISH_distribution.pdf"), 
       plot=main, 
       width = 14,
       height = 12)

# ---------------------------------------------------------------------------------
