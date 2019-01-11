# This is the main R file to generate the results related our manuscript
# figure 2 and the supplementary figures related to figure 2.
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


# the following script is to generate the simulation data. Here we use 
# HPC to generate the simulation which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    generate_save_data(dropout_index, seed_value)
    
  }
}

# The reader could skip the above script and we have deposit the genereate 
# data in the current folder /simulation_data/

# the following codes are used to do the imputation using different methods

# the following script is to impute data. Here we use 
# HPC to impute the data using drimpute which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_drimpute(dropout_index, seed_value)
    
  }
  
}

# the following script is to impute data. Here we use 
# HPC to impute the data using scimpute which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_scimpute(dropout_index, seed_value)
    
  }
  
}

# the following script is to impute data. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    run_scrabble(dropout_index, seed_value)
    
  }
}

# the following script is to calculate the errors. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:100)){
    
    result <- calculate_error_splatter(drop_index, seed_value)
    
    dir.create(file.path('error_data'), showWarnings = FALSE)
    
    saveRDS(result,
            file = paste0("error_data/error_",drop_index,"_",seed_value,".rds")
    )
    
  }
  
}

# Gather the errors
# -----------------
error_list <- list()

error_cell_list <- list()

error_gene_list <- list()

# gather the error from the data from different dropout rates
for(i in c(1:3)){
  error_matrix <- c()
  
  error_cell_matrix <- c()
  
  error_gene_matrix <- c()
  
  for(j in c(1:100)){
    
    tmp <- readRDS(file = paste0("error_data/error_",i,"_",j,".rds"))
    
    error_matrix <- cbind(error_matrix,as.matrix(tmp$error))
    
    error_cell_matrix <- cbind(error_cell_matrix,as.matrix(tmp$error_cell))
    
    error_gene_matrix <- cbind(error_gene_matrix,as.matrix(tmp$error_gene))
    
  }
  
  error_list[[i]] <- error_matrix
  
  error_cell_list[[i]] <- error_cell_matrix
  
  error_gene_list[[i]] <- error_gene_matrix
  
}

# save the errors
saveRDS(error_list,file = "error_all.rds")

saveRDS(error_cell_list,file = "error_all_cell.rds")

saveRDS(error_gene_list,file = "error_all_gene.rds")
# ------

# Plot the boxplots in figure 2 and the supplementary figures related to Figure 2
# ----------------------------------------------------------------------------
# load the error data
error_list <- readRDS(file = "error_all.rds")

error_cell_list <- readRDS(file = "error_all_cell.rds")

error_gene_list <- readRDS(file = "error_all_gene.rds")


# define the list for boxplot comparisons
p <- list()

# Dropout rate: 71%
p[[1]] <- plot_comparison(error_list[[1]], "Error", 1400, 140)

p[[2]] <- plot_comparison(error_cell_list[[1]], "Correlation", 1, 0.1)

p[[3]] <- plot_comparison(error_gene_list[[1]], "Correlation", 0.4, 0.04)

# Dropout rate: 83%
p[[4]] <- plot_comparison(error_list[[2]], "Error", 1800, 180)

p[[5]] <- plot_comparison(error_cell_list[[2]], "Correlation", 1, 0.1)

p[[6]] <- plot_comparison(error_gene_list[[2]], "Correlation", 0.3, 0.03)

# Dropout rate:87%
p[[7]] <- plot_comparison(error_list[[3]], "Error",1800,180)

p[[8]] <- plot_comparison(error_cell_list[[3]], "Correlation", 1, 0.1)

p[[9]] <- plot_comparison(error_gene_list[[3]], "Correlation", 0.2, 0.02)

# save the PDF files
main <- grid.arrange(grobs = p,ncol = 3)
ggsave(filename="Figure_error.pdf", 
       plot = main, 
       width = 18, 
       height = 12)
# ----------------------------------------------------------------------------


# plot the mean-variance figures in Figure 2 and the supplementary figures realted to Figure 2
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_meansd(drop_index, seed_value)
  
  ggsave(filename=paste0("Figures_mean_variance_",drop_index,".pdf"),
         plot = p,
         width = 24,
         height = 3)
}
# ----------------------------------------------------------------------------


# plot the mean-variance figures in Figure 2 and the supplementary figures realted to Figure 2
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_comparison_tsne(drop_index,
                            seed_value,
                            50,100)
  
  ggsave(filename=paste0("Figures_tsne_",drop_index,".pdf"),
         plot = p,
         width = 24,
         height = 3)
}
# -------------------------------------------------------------------------

# plot the mean-variance figures in Figure 2 and the supplementary figures realted to Figure 2
# ----------------------------------------------------------------------------
for(drop_index in c(1:3)){
  
  seed_value <- 10
  
  p <- plot_ma(drop_index, seed_value)
  
  ggsave(filename=paste0("Figure_ma_",drop_index,".pdf"), 
         plot = p, 
         width = 24, 
         height = 3)
}
# -------------------------------------------------------------------------


