# This is the main R file to generate the results related our manuscript
# figure 3 and the supplementary figures related to figure 3.
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

# Load the data
# ----------------------------------------------------------------------------

# load the bulk RNAseq data, sample informaiton, and the gene name information
data_bulk <- read.table(file = "data/data_bulk.txt", header = TRUE)

sample_name_df <- read.table(file = "data/sample_info.csv", header = FALSE)
cell_type <- unique(as.character(sample_name_df$V1))
sample_name <- as.character(sample_name_df$V1)

gene_name <- read.table(file = "data/mouse_gene_id.txt", sep = ",")

# filter out the genes with low expression 
index <- rowSums(name_file > 1) > 1
data <- log10(name_file[index,] + 1)
gene_name <- gene_name[index,]

# To reduce the computation cost, we sample 5000 genes from those all genes
data1 <- as.matrix(data[sample(1:dim(data)[1],5000),])
# ----------------------------------------------------------------------------


# Prepare the true single cell data
# ----------------------------------------------------------------------------
data_noise <- prepare_data(data1, cell_type, sample_name)

# save all the data 
save(data,data1, cell_type, sample_name, data_noise, file = "data_raw.RData")
# ----------------------------------------------------------------------------

# the following script is to generate the simulation data. Here we use 
# HPC to generate the simulation which could reduce the running time
# ----------------------------------------------------------------------------
for(dropout_index in c(1:4) ){
  for(seed_value in c(1:100)){
    
    run_generate_data(dropout_index, seed_value)
    
  }
}
# ----------------------------------------------------------------------------

# the following script is to run DrImpute. Here we use 
# HPC to generate the simulation which could reduce the running time
# ----------------------------------------------------------------------------
for(dropout_index in c(1:4) ){
  for(seed_value in c(1:100)){
    
    run_drimpute(dropout_index, rand_seed)
    
  }
}
# ----------------------------------------------------------------------------

# the following script is to run scImpute. Here we use 
# HPC to generate the simulation which could reduce the running time
# ----------------------------------------------------------------------------
for(dropout_index in c(1:4) ){
  for(seed_value in c(1:100)){
    
    run_scimpute(dropout_index, rand_seed)
    
  }
}
# ----------------------------------------------------------------------------

# the following script is to run SCRABBLE. Here we use 
# HPC to generate the simulation which could reduce the running time
# ----------------------------------------------------------------------------
for(dropout_index in c(1:4) ){
  for(seed_value in c(1:100)){
    
    run_scrabble(dropout_index, rand_seed)
    
  }
}
# ----------------------------------------------------------------------------

# the following script is to calculate the error. Here we use 
# HPC to generate the simulation which could reduce the running time
# ----------------------------------------------------------------------------
for(dropout_index in c(1:4) ){
  for(seed_value in c(1:100)){
    
    result <- run_error(dropout_index, rand_seed)
    
    dir.create(file.path('error_data'), showWarnings = FALSE)
    
    saveRDS(result,
            file = paste0("error_data/error_",
                          drop_index,"_",
                          seed_value,".rds")
    )
    
  }
}
# ----------------------------------------------------------------------------

# the following script is to assemble the errors. Here we use 
# HPC to generate the simulation which could reduce the running time
# ----------------------------------------------------------------------------
error_list <- list()
error_cell_list <- list()
error_gene_list <- list()

for(i in c(1:4)){
  
  error_matrix <- c()
  error_cell_matrix <- c()
  error_gene_matrix <- c()
  
  for(j in c(1:100)){
    
    # load the data
    tmp <- readRDS(file = paste0("error_data/error_",i,"_",j,".rds"))
    
    error_matrix <- cbind(error_matrix,as.matrix(tmp$error))
    
    error_cell_matrix <- cbind(error_cell_matrix,as.matrix(tmp$error_cell))
    
    error_gene_matrix <- cbind(error_gene_matrix,as.matrix(tmp$error_gene))
    
  }
  
  error_list[[i]] <- error_matrix
  error_cell_list[[i]] <- error_cell_matrix
  error_gene_list[[i]] <- error_gene_matrix
  
}

saveRDS(error_list,file = "error_all.rds")
saveRDS(error_cell_list,file = "error_all_cell.rds")
saveRDS(error_gene_list,file = "error_all_gene.rds")
# ----------------------------------------------------------------------------

# Plot the boxplot with pvalues for the errors
# ----------------------------------------------------------------------------
# load the data
error_list <- readRDS(file = "error_all.rds")
error_cell_list <- readRDS(file = "error_all_cell.rds")
error_gene_list <- readRDS(file = "error_all_gene.rds")

# define the graph list
p <- list()

p[[1]] <- plot_comparison(error_list[[1]], "Error", 1200, 80)
p[[2]] <- plot_comparison(error_cell_list[[1]], "Correlation", 1,0.1)
p[[3]] <- plot_comparison(error_gene_list[[1]], "Correlation",1,0.1)
p[[4]] <- plot_comparison(error_list[[2]], "Error",1200,80)
p[[5]] <- plot_comparison(error_cell_list[[2]], "Correlation",1,0.1)
p[[6]] <- plot_comparison(error_gene_list[[2]], "Correlation",1,0.1)
p[[7]] <- plot_comparison(error_list[[3]], "Error",1200,80)
p[[8]] <- plot_comparison(error_cell_list[[3]], "Correlation",1,0.1)
p[[9]] <- plot_comparison(error_gene_list[[3]], "Correlation",1,0.1)
p[[10]] <- plot_comparison(error_list[[4]], "Error",1200,80)
p[[11]] <- plot_comparison(error_cell_list[[4]], "Correlation",1,0.1)
p[[12]] <- plot_comparison(error_gene_list[[4]], "Correlation",1,0.1)

# save the figures
main <- grid.arrange(grobs = p,ncol = 3)
ggsave(filename="Figure_HF_error.pdf",
       plot = main,
       width = 18,
       height = 15)
# ----------------------------------------------------------------------------


# Plot tSNE
# ----------------------------------------------------------------------------
# define the random seed
rand_seed <- 2

# define the plot list
p <- list()

# plot the tSNE
for(dropout_index in c(1:4)){
  
  p[[dropout_index]] <- plot_tsne_HF(dropout_index, rand_seed , k = 4)
  
}

# Assemble all figures in one graph
main <- grid.arrange(grobs = p,ncol = 1)

# Save the plots
ggsave(
  filename = paste0("Figure_HF_tsne.pdf"),
  plot = main,
  width = 30,
  height = 16
)
# ----------------------------------------------------------------------------


# Plot MA
# ----------------------------------------------------------------------------
# define the random seed
rand_seed <- 2

# define the plot list
p <- list()

# plot the tSNE
for(dropout_index in c(1:4)){
  
  p[[dropout_index]] <- plot_ma_HF(dropout_index, rand_seed , k = 4)
  
}

# Assemble all figures in one graph
main <- grid.arrange(grobs = p,ncol = 1)

# Save the plots
ggsave(
  filename = paste0("Figure_HF_ma.pdf"),
  plot = main,
  width = 30,
  height = 16
)
# ----------------------------------------------------------------------------

# Plot Mean_Variance
# ----------------------------------------------------------------------------
# define the random seed
rand_seed <- 2

# define the plot list
p <- list()

# plot the tSNE
for(dropout_index in c(1:4)){
  
  p[[dropout_index]] <- plot_meansd(dropout_index, rand_seed , k = 4)
  
}

# Assemble all figures in one graph
main <- grid.arrange(grobs = p,ncol = 1)

# Save the plots
ggsave(
  filename = paste0("Figure_HF_meansd.pdf"),
  plot = main,
  width = 30,
  height = 16
)
# ----------------------------------------------------------------------------


# Plot data
# ----------------------------------------------------------------------------
# define the random seed
rand_seed <- 2

# define the plot list
p <- list()

# plot the tSNE
for(dropout_index in c(1:4)){
  
  p[[dropout_index]] <- plot_data(dropout_index, rand_seed)
  
}

# Assemble all figures in one graph
main <- grid.arrange(grobs = p,ncol = 1)

# Save the plots
ggsave(
  filename = paste0("Figure_HF_meansd.pdf"),
  plot = main,
  width = 30,
  height = 16
)
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# define the random seed
rand_seed <- 2

# define the plot list
p <- list()

# plot the tSNE
for(dropout_index in c(1:4)){
  
  p <- plot_cor_HF(dropout_index, rand_seed)
  
  ggsave(
    filename = paste0("HF_Cor_Gene_",lambda_index,".pdf"),
    plot = p[[1]],
    width = 35,
    height = 4
  )
  
  ggsave(
    filename = paste0("HF_Cor_Cell_",lambda_index,".pdf"),
    plot = p[[2]],
    width = 35,
    height = 4
  )
  
}
# ----------------------------------------------------------------------------









