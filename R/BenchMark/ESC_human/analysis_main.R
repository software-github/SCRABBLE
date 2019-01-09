# This is the main R file to generate the results related our manuscript
# figure 6 and the supplementary figures related to figure 6.
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

# load the data
# ---------------------------------------------------------------------------------
data_sc <- fread(file = "data_all/data_sc_TPM.csv", 
                 sep=',')

data_bulk <- fread(file = "data_all/data_bulk_TPM.csv",
                   sep = ',')

sample_sc_information <- fread(file = "data_all/single_sample_information.csv",
                               sep='\t')

sample_bulk_information <- fread(file  = "data_all/bulk_sample_information.csv",sep='\t', 
                                 header = FALSE)

gene_name <- fread(file = "data_all/gene_name.csv",
                   sep=',')


# ---------------------------------------------------------------------------------

# Preprocess the data
# ---------------------------------------------------------------------------------
celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

data_sc <- as.data.frame(data_sc)

data_bulk <- as.data.frame(data_bulk)

for(i in c(1:7)){
  
  index <- sample_sc_information$V1 == celltype[i]
  
  print(paste0(celltype[i],': ',sum(index)))
  
  data_sc_tmp <- data_sc[,which(index)]
  
  index <- sample_bulk_information$V1 == celltype[i]
  
  data_bulk_tmp <- data_bulk[,index]
  
  tmp_sc <- rowSums(data_sc_tmp > 0) == dim(data_sc_tmp)[2]
  
  mean_sc <- rowMeans(data_sc_tmp)
  
  mean_bulk <- rowMeans(data_bulk_tmp)
  
  x <- mean_sc[tmp_sc]
  
  y <- mean_bulk[tmp_sc]
  
  c1 <- lm(x ~ y + 0)
  
  data_bulk_tmp <- data_bulk_tmp*c1$coefficients
  
  fwrite(data_sc_tmp, 
         file = paste0("data_all/data_sc_", celltype[i],".csv"),
         sep = ',')
  
  fwrite(data_bulk_tmp,
         file = paste0("data_all/data_bulk_", celltype[i],".csv"),
         sep = ',')
} 

# ---------------------------------------------------------------------------------
# Run DrImpute 
# ---------------------------------------------------------------------------------

dir.create(file.path("/imputation_drimpute_data/"))

celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

for(i in c(1:7)){
  
  data_sc <- as.matrix(fread(paste0(cwd,"/data_all/data_sc_",celltype[i],".csv")))
  
  data_bulk <- as.matrix(fread(paste0(cwd,"/data_all/data_bulk_",celltype[i],".csv")))
  
  load(file = paste0("data_all/IPA",
                     "_",
                     celltype[i],
                     "_index.RData")
  )
  
  data_sc1 <- data_sc
  
  data_sc <- data_sc[index_sc,]
  
  extdata <- DrImpute(log10(data_sc + 1))
  
  extdata <- 10^extdata - 1
  
  data_sc1[index,] <- extdata
  
  saveRDS(data_sc1, 
          file = paste0("/imputation_drimpute_data/data_",celltype[i],"_drimpute_imputation.rds"))
  
}

# ---------------------------------------------------------------------------------
# Run scImpute 
# ---------------------------------------------------------------------------------

celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

for(i in c(1:7)){
  
  data_sc <- as.matrix(fread(paste0("/data_all/data_sc_",celltype[i],".csv")))
  
  data_bulk <- as.matrix(fread(paste0("/data_all/data_bulk_",celltype[i],".csv")))
  
  load(file = paste0("data_all/IPA",
                     "_",
                     celltype[i],
                     "_index.RData")
  )
  
  data_sc1 <- data_sc
  
  data_sc <- data_sc[index_sc,]
  
  dir.create(file.path("/imputation_scimpute_data/"), showWarnings = FALSE)
  
  # impute the data using scImpute
  data_dropout <- data_sc
  
  write.table(data_dropout,paste0("/imputation_scimpute_data/sdropout_scimpute_",i,".csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  # run scImpute to obtain the missing data
  scimpute(
    paste0("/imputation_scimpute_data/dropout_scimpute_",i,".csv"),
    infile = "csv",           
    outfile = "csv",          
    out_dir = paste0(cwd1, "scImpute_", i, "_"),
    drop_thre = 0.5,          
    Kcluster = 2,
    ncores = 1)     
  # 
  # # clean the data
  data_dropout <- read.table( file = paste0("/imputation_scimpute_data/scImpute_",
                                            i,
                                            "_scimpute_count.csv") ,
                              header = TRUE, sep=",")
  
  data_dropout$X <- NULL
  
  data_sc1[index_sc,] <- as.matrix(data_dropout)
  
  saveRDS(data_sc1, 
          file = paste0("/imputation_scimpute_data/data_",celltype[i],"_scimpute_imputation.rds"))
  
}

# ---------------------------------------------------------------------------------
# Run MAGIC 
# ---------------------------------------------------------------------------------

# import sys
# import os
# import magic 
# import pandas as pd
# 
# cwd = os.getwd()
# print cwd
# if not os.path.exists(cwd+"/imputation_magic_data"):
#   os.makedirs(cwd+"/imputation_magic_data")
# 
# celltype = ["H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1"]
# 
# 
# threshold = 1e-10
# 
# for i in range(7):
#
#     X = pd.read_csv(cwd + "data_all/data_sc_"+celltype[i]+".csv",sep = ',')
# 
#     Y = pd.read_csv(cwd + "data_all/data_bulk_"+celltype[i]+".csv",sep = ',')
# 
#     X1 = X
# 
#     X = X.loc[X.var(axis=1) > threshold]
# 
#     magic_operator = magic.MAGIC(n_pca = 50)
# 
#     X_magic = magic_operator.fit_transform(X.T)
# 
#     out_magic = X_magic.T
# 
#     X1.loc[X1.var(axis=1) > threshold] = out_magic
# 
#     X1.to_csv(cwd+"/imputation_magic_data/data_magic_"+celltype[i]+".csv", sep = ',', header= None)

# ---------------------------------------------------------------------------------
# Run SCRABBLE
# ---------------------------------------------------------------------------------

dir.create(file.path("/imputation_scrabble_data/"))

celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

for(i in c(1:7)){
  
  data_sc <- as.matrix(fread(paste0("/data_all/data_sc_",celltype[i],".csv")))
  
  data_bulk <- as.matrix(fread(paste0("/data_all/data_bulk_",celltype[i],".csv")))
  
  
  load(file = paste0("data_all/IPA",
                     "_",
                     celltype[i],
                     "_index.RData")
  )
  
  data_sc1 <- data_sc
  
  data_sc <- data_sc[index_sc,]
  
  data_bulk <- rowMeans(data_bulk[index_sc,])
  
  data <- list()
  
  data[[1]] <- data_sc
  
  data[[2]] <- data_bulk
  
  parameter <- c(1,1e-5,1e-4)
  
  print(parameter)
  
  result <- scrabble(data,
                     parameter = parameter, 
                     nIter = 30,
                     error_out_threshold = 1e-14, 
                     nIter_inner = 30,
                     error_inner_threshold = 1e-14)
  
  data_sc1[index_sc,] <- result
  
  saveRDS(data_sc1, file = paste0("/imputation_scrabble_data/data_",
                                  celltype[i],"_scrabble_imputation.rds"))
  
}


# Generate the index of the pathway databases

# ---------------------------------------------------------------------------------
# Run IPA
# ---------------------------------------------------------------------------------

# load the gene
data_gene <- fread(file = "data_all/gene_ESC.csv", 
                   header = FALSE)

data_sc <- as.matrix(fread(file = paste0("data_all/data_sc_",
                                         "DEC",
                                         ".csv")))


var0 <- MatVar(data_sc,1)

index_sc <- var0 > 1e-10

data_gene <- data_gene[index_sc,]


n_gene <- dim(data_gene)[1]


# IPA
N <- 186
index <- list()
for(i in c(1:N)){
  
  tmp <- fread(file = paste0("data_all/IPA/IPA_gene_",i,".csv"), 
               header = FALSE)
  
  tmp <- match(tmp$V1, data_gene$V1)
  
  tmp <- tmp[!is.na(tmp)]
  
  index[[i]] <- tmp
}

for(i in c(1:N)){
  
  tmp <- c()
  
  k <- length(index[[i]])
  
  for(j in c(1:100)){
    set.seed(j)
    tmp <- rbind(tmp, sample(1:n_gene,k))
  }
  
  index[[i]] <- rbind(index[[i]],tmp)
  
}

saveRDS(index, file = "data_all/IPA_index.rds")

# ---------------------------------------------------------------------------------
# RunKEGG
# ---------------------------------------------------------------------------------

N <- 186
index <- list()
for(i in c(1:N)){
  
  tmp <- fread(file = paste0("data_all/KEGG/KEGG_gene_",i,".csv"), 
               header = FALSE)
  
  tmp <- match(tmp$V1, data_gene$V1)
  
  tmp <- tmp[!is.na(tmp)]
  
  index[[i]] <- tmp
}

for(i in c(1:N)){
  
  tmp <- c()
  
  k <- length(index[[i]])
  
  for(j in c(1:100)){
    set.seed(j)
    tmp <- rbind(tmp, sample(1:n_gene,k))
  }
  
  index[[i]] <- rbind(index[[i]],tmp)
  
}

saveRDS(index, file = "data_all/KEGG_index.rds")

# ---------------------------------------------------------------------------------
# Run REACTOME
# ---------------------------------------------------------------------------------
# 
N <- 674
index <- list()
for(i in c(1:N)){
  
  tmp <- fread(file = paste0("data_all/REACTOME/REACTOME_gene_",i,".csv"), 
               header = FALSE)
  
  tmp <- match(tmp$V1, data_gene$V1)
  
  tmp <- tmp[!is.na(tmp)]
  
  index[[i]] <- tmp
}

for(i in c(1:N)){
  
  tmp <- c()
  
  k <- length(index[[i]])
  
  for(j in c(1:100)){
    set.seed(j)
    tmp <- rbind(tmp, sample(1:n_gene,k))
  }
  
  index[[i]] <- rbind(index[[i]],tmp)
  
}

saveRDS(index, file = "data_all/REACTOME_index.rds")

# ---------------------------------------------------------------------------------
#  Generate the index 
# ---------------------------------------------------------------------------------

celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

pathways <- c("IPA", "KEGG", "REACTOME")

for(i in c(1:7)){
  for(j in c(1:3)){
    
    generate_index(celltype[i], pathways[j])
    
  }
}

# ---------------------------------------------------------------------------------
#  Generate the ratio
# ---------------------------------------------------------------------------------


celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

pathways <- c("IPA", "KEGG", "REACTOME")

method_name <- c("dropout","drimpute","scimpute","magic","scrabble")

for(j in c(1:3)){

  for(i in c(1:7)){
    
    for(k in c(1:5)){
      
      values <- calculate_ratio(celltype[i], pathways[j], method_name[k])
      
      saveRDS(values, file = paste0(cwd,"/data_all/data_",
                                    celltype[i], "_",pathways[j],
                                    "_", method_name[k],".rds"))
    }
    
  }
  
}

# ---------------------------------------------------------------------------------
#  Plot the ratio
# ---------------------------------------------------------------------------------

celltype <- c("H9", "DEC", "EC",  "HFF", "NPC", "TB", "H1")

pathways <- c("IPA", "KEGG", "REACTOME")

method_name <- c("dropout","drimpute","scimpute","magic","scrabble")

for(j in c(1:3)){
  
  pl <- list()
  
  for(i in c(1:7)){

    for(k in c(1:5)){
      
      values <- readRDS(file = paste0(cwd,"/data_all/data_",
                                      celltype[i], "_",pathways[j],
                                      "_", method_name[k],".rds"))
      
      ratio <- 100*((values[,1]/values[,2]) - 1)
      
      tmp <- rbind(tmp, ratio)
      
    }
    
    pl[[i]] <- plot_pathway_ratio(tmp, celltype[i], pathways[j])
    
  }
  
  main <- grid.arrange(grobs = pl,ncol = 7, top = "main")
  
  ggsave(
    filename = paste0("Figure_pathway_",pathways[j],".pdf"),
    plot = main,
    width = 30,
    height = 5
  )
  
}