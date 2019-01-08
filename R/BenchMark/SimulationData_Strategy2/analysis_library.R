# generate the dropout data
run_generate_data <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the data
  load("data_raw.RData")
  
  # the dropout parameter
  lambda = c(0.3,0.4,0.6,0.8)
  
  # select the index of cells with 800 cells from 2000 cells
  index_matrix = t(matrix(1:2000,ncol = 20))
  
  set.seed(rand_seed)
  
  index_sample = as.vector(t(index_matrix[sample(1:20,8),]))
  
  # define the true data
  data_true = data_noise[,index_sample]
  
  # define the dropout 
  p = matrix(rep(exp(-lambda[dropout_index]*rowMeans(data_true)^2),
                 dim(data_true)[2]), ncol = dim(data_true)[2])
  
  set.seed(rand_seed)
  
  dropout_s = matrix(runif(length(p)), nrow = dim(p)[1])
  
  dropmatrix = dropout_s > p 
  
  # define the raw data
  data_raw = data_true*dropmatrix
  
  # define the bulk data
  data_bulk = rowMeans(data_true)
  
  # save the raw data
  fwrite(as.data.frame(data_raw), 
         file = paste0(cwd,"/data_all/data_raw_",dropout_index,"_",rand_seed,".csv"))
  
  # save the true data
  fwrite(as.data.frame(data_true), 
         file = paste0(cwd,"/data_all/data_true_",dropout_index,"_",rand_seed,".csv"))
  
  # save the bulk data
  fwrite(as.data.frame(data_bulk), 
         file = paste0(cwd,"/data_all/data_bulk_",dropout_index,"_",rand_seed,".csv"))
  
}

# prepare the true single cell data
prepare_data <- function(data1, cell_type, sample_name){
  
  # Parameter in the function
  # data1: the bulk RNAseq data with 5000 genes
  # cell_type: the cell type names
  # sample_name: all the sample information
  
  
  # expand the data
  data_generate = c()
  
  for (i in c(1:length(cell_type))){
    
    data_generate = cbind(data_generate,rowMeans(data1[,sample_name == cell_type[i]]))
    
  }
  
  # add the noise
  k_seed = 1
  
  data_noise = data_generate[,rep(1:20,each = 100)]
  
  # add the noise to the data
  for (i in c(1:length(cell_type))){
    
    var = apply(data1[,sample_name == cell_type[i]],1,sd)
    
    for (j in c(1:length(var))){
      
      set.seed(k_seed)
      
      k_seed = k_seed + 1
      
      data_noise[j,(100*(i - 1) + 1):(100*(i - 1) + 100)] = data_noise[j,(100*(i - 1) + 1):(100*(i - 1) + 100)] + 
        t(rnorm(100,0,5*var[j]))
    }
  }
  
  # avoid the negative values and the too large expression values
  data_noise[data_noise < 0] = 0
  data_noise[data_noise > 5] = 5
  
  return(data_noise)
  
}

# run DrImpute 
run_drimpute <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # create the folder
  dir.create(file.path("/imputation_drimpute_data/"))
  
  # load the data
  data_sc = as.matrix(fread(paste0(cwd,"/data_all/data_raw_",dropout_index,"_",rand_seed,".csv")))
  
  # run the imputation
  extdata = DrImpute(data_sc)
  
  # save the data
  saveRDS(extdata, file = paste0(cwd1,"data_",dropout_index,"_",rand_seed,"_drimpute_imputation.rds"))
  
}

# run scImpute 
run_scimpute <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # create the folder
  dir.create(file.path("/imputation_scimpute_data/"), showWarnings = FALSE)
  dir.create(file.path("/temp_scimpute_data/"), showWarnings = FALSE)
  
  # load the data
  data_sc = as.matrix(fread(paste0("/data_all/data_raw_",
                                   dropout_index,"_",
                                   rand_seed,".csv")))
  
  # impute the data using scImpute
  write.table(data_sc,
              paste0("temp_scimpute_data/dropout_scimpute_",
                     dropout_index,"_",
                     rand_seed,".csv"),
              sep=',',
              row.names = TRUE,
              col.names = TRUE
  )
  
  # run scImpute to obtain the missing data
  scimpute(
    paste0("temp_scimpute_data/dropout_scimpute_",dropout_index,"_",rand_seed,".csv"),
    infile = "csv",          
    outfile = "csv",          
    out_dir = paste0("temp_scimpute_data/scImpute_", dropout_index, "_", rand_seed,"_"),   
    drop_thre = 0.5,          
    Kcluster = 2,
    ncores = 2)         

  # clean the data
  data_scimpute = read.table( file = paste0("temp_scimpute_data/scImpute_",
                                             dropout_index, "_",
                                             rand_seed, "_scimpute_count.csv"),
                              header = TRUE, sep=",")
  data_scimpute$X = NULL
  
  # save the data
  saveRDS(data_scimpute, file = paste0("/imputation_scimpute_data/data_",
                                       dropout_index, "_",
                                       rand_seed, "_scimpute_imputation.rds"))

}

# run MAGIC
# The following is the python script used in the analysis
# ----------------------------------------------------------
# import sys
# import os
# import magic 
# import pandas as pd
# 
# cwd = os.getcwd()
# 
# if not os.path.exists(cwd+"/magic_data"):
#   os.makedirs(cwd+"/magic_data")
# 
# X =pd.read_csv(cwd + "/data_all/data_raw_"+str(dropout_value)+"_"+str(seed_value)+".csv",sep = ',')
# magic_operator = magic.MAGIC()
# X_magic = magic_operator.fit_transform(X.T)
# out_magic = X_magic.T
# out_magic.to_csv(cwd+"/magic_data/magic_"+str(dropout_value)+"_"+str(seed_value)+".csv", sep = ',', header= None)
# ----------------------------------------------------------

# run SCRABBLE
run_scrabble <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # create the folder
  dir.create(file.path("/imputation_scrabble_data/"), showWarnings = FALSE)
  
  # load the data
  data_sc = as.matrix(fread(paste0("/data_all/data_raw_",
                                    dropout_index,"_",
                                    rand_seed,".csv")))
  
  data_bulk = as.matrix(fread(paste0("/data_all/data_bulk_",
                                      dropout_index,"_",
                                      rand_seed,".csv")))
  
  # Prepare the data
  data = list()
  data[[1]] = data_sc
  data[[2]] = data_bulk
  
  # parameter setting
  parameter = c(1,1e-5,1e-2)
  
  # run scrabble
  result = scrabble(data,
                    parameter = parameter, 
                    nIter = 30,
                    error_out_threshold = 1e-7, 
                    nIter_inner = 30,
                    error_inner_threshold = 1e-5)
  
  # save the data
  saveRDS(result, file = paste0(cwd1,"data_",i,"_",j,"_",k,"_scrabble_imputation.rds"))
  
}

# get the data, true data, raw data, and the imputed data
get_data_HF <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load the true data
  data_true = as.matrix(fread(paste0("/data_all/data_true_",
                                     dropout_index,"_",
                                     rand_seed,".csv")))
  
  # load the dropout data
  data_raw = as.matrix(fread(paste0("/data_all/data_raw_",
                                    dropout_index,"_",
                                    rand_seed,".csv")))
  
  # load the imputed data by DrImpute
  data_drimpute = readRDS(file = paste0("/imputation_drimpute_data/data_",
                                        dropout_index,"_",
                                        rand_seed,"_drimpute_imputation.rds"))
  
  # load the imputed data by scImpute
  data_scimpute = readRDS(file = paste0(cwd,"/imputation_scimpute_data/data_",
                                        dropout_index,"_",
                                        rand_seed,"_scimpute_imputation.rds"))
  
  # load the imputed data by MAGIC
  data_magic = as.matrix(fread(paste0("/magic_data/magic_",
                                      dropout_index,"_",
                                      rand_seed,".csv")))
  data_magic = data_magic[,-1]
  
  # load the imputed data by SCRABBLE
  data_scrabble = readRDS(file = paste0("/imputation_scrabble_data/data_",
                                        dropout_index,"_",
                                        rand_seed,"_scrabble_imputation.rds"))
  

  # define the data as a list
  data = list()
  
  data$data_true = as.matrix(data_true)
  
  data$data_raw = as.matrix(data_raw)
  
  data$data_drimpute = as.matrix(data_drimpute)
  
  data$data_scimpute = as.matrix(data_scimpute)
  
  data$data_magic = as.matrix(data_magic)
  
  data$data_scrabble = as.matrix(data_scrabble)
  
  return(data)
  
}

# get the similarity between two datasets
calculate_similarity <- function(data1,data2){
  
  d = cor(c(data1[lower.tri(data1)]),c(data2[lower.tri(data2)]))
  
  return(d)
  
}

# calculate the error
calculate_error <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load all the data
  data = get_data_HF(dropout_index, rand_seed)
  
  data_true  = data$data_true
  data_dropout = data$data_raw
  data_drimpute = data$data_drimpute
  data_scimpute = data$data_scimpute
  data_magic = data$data_magic
  data_scrabble = data$data_scrabble
  
  error = matrix(0, nrow = 6, ncol = 1)
  error[1] = norm(data_dropout - data_true, type = "2")
  error[2] = norm(data_drimpute - data_true, type = "2")
  error[3] = norm(data_scimpute - data_true, type = "2")
  error[4] = norm(data_magic - data_true, type = "2")
  error[5] = norm(data_scrabble - data_true, type = "2")
  error[6] = 1 - nnzero(data_dropout)/length(data_dropout)
  
  # gene-gene correlation
  # true data
  data_true_gene = cor(t(data_true), method = "pearson")
  data_true_gene[is.na(data_true_gene)] = 0
  
  # dropout data
  data_dropout_gene = cor(t(data_dropout), method = "pearson")
  data_dropout_gene[is.na(data_dropout_gene)] = 0
  
  # DrImpute data
  data_drimpute_gene = cor(t(data_drimpute), method = "pearson")
  data_drimpute_gene[is.na(data_drimpute_gene)] = 0
  
  # scImpute data
  data_scimpute_gene = cor(t(data_scimpute), method = "pearson")
  data_scimpute_gene[is.na(data_scimpute_gene)] = 0
  
  # MAGIC data
  data_magic_gene = cor(t(data_magic), method = "pearson")
  data_magic_gene[is.na(data_magic_gene)] = 0
  
  # SCRABBLE data
  data_scrabble_gene = cor(t(data_scrabble), method = "pearson")
  data_scrabble_gene[is.na(data_scrabble_gene)] = 0
  
  
  error_gene = matrix(0, nrow = 6, ncol = 1)
  error_gene[1] = calculate_similarity(data_true_gene, data_dropout_gene) 
  error_gene[2] = calculate_similarity(data_true_gene, data_drimpute_gene)
  error_gene[3] = calculate_similarity(data_true_gene, data_scimpute_gene)
  error_gene[4] = calculate_similarity(data_true_gene, data_magic_gene)
  error_gene[5] = calculate_similarity(data_true_gene, data_scrabble_gene)
  error_gene[6] = 1 - nnzero(data_dropout)/length(data_dropout)
  
  
  # cell-cell correlation
  
  # true data
  data_true_cell = cor(data_true, method = "pearson")
  data_true_cell[is.na(data_true_cell)] = 0
  
  # dropout data
  data_dropout_cell = cor(data_dropout, method = "pearson")
  data_dropout_cell[is.na(data_dropout_cell)] = 0
  
  # DrImpute data
  data_drimpute_cell = cor(data_drimpute, method = "pearson")
  data_drimpute_cell[is.na(data_drimpute_cell)] = 0
  
  # scImpute data
  data_scimpute_cell = cor(data_scimpute, method = "pearson")
  data_scimpute_cell[is.na(data_scimpute_cell)] = 0
  
  # MAGIC data
  data_magic_cell = cor(data_magic, method = "pearson")
  data_magic_cell[is.na(data_magic_cell)] = 0
  
  # SCRABBLE
  data_scrabble_cell = cor(data_scrabble, method = "pearson")
  data_scrabble_cell[is.na(data_scrabble_cell)] = 0
  
  error_cell = matrix(0, nrow = 6, ncol = 1)
  error_cell[1] = calculate_similarity(data_true_cell, data_dropout_cell)
  error_cell[2] = calculate_similarity(data_true_cell, data_drimpute_cell)
  error_cell[3] = calculate_similarity(data_true_cell, data_scimpute_cell)
  error_cell[4] = calculate_similarity(data_true_cell, data_magic_cell)
  error_cell[5] = calculate_similarity(data_true_cell, data_scrabble_cell)
  error_cell[6] = 1 - nnzero(data_dropout)/length(data_dropout)
  
  # define the error as a list
  result <- list()
  result$error <- error
  result$error_cell <- error_cell
  result$error_gene <- error_gene
  
  return(result)
  
}

# plot the comparison of the errors
plot_comparison <- function(data, ylables, ylim_value = 100, h_ylim = 10){
  
  # this function is used to plot the boxplot with the p-values
  # Parameter in the function
  # data : data is a matrix with 6 columns, the first five colums are the errors
  # the last column is the dropout rates
  # ylabels: the y labels shown on the graph
  # ylim_value : the smallest y value where the pvalue is shown
  # h_ylim: the difference between two pvalues shown in the graph. the best ones
  # is the 10%*ylim_vlaue
  
  
  dataV0 = data[c(1:5),]
  
  dataV1 = data.frame(as.vector(t(dataV0)))
  
  # calculate the dropout rate
  dropout_rate = round(mean(data[6,])*100)
  
  N = dim(dataV1)[1]                                                      
  
  # the number of data using for plotting the boxplot
  # dataV1 is built with two columns: y values and group labels
  dataV1$group = rep(c(1:5), each = N/5)
  
  colnames(dataV1) = c('y','group')
  
  # define the comparison lists
  my_comparisons = list( c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))
  
  # compare the errors using t-test
  pval = compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "5", paired = TRUE)

  # plot the boxplot with pvalues for the comparisons
  pp <- ggboxplot(dataV1, x = "group", y = "y", fill = "group",
                  palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00")) +
    stat_boxplot(geom = "errorbar", width = 0.3, outlier.size = NA, outlier.shape = NA, outlier.colour = NA) +
    ylim(c(0,ylim_value + 3.5*h_ylim)) + 
    theme_bw() +
    geom_signif(comparisons = my_comparisons, 
                annotations = formatC(pval$p, format = "e", digits = 2),
                tip_length = 0.01,
                y_position = c(ylim_value + 3*h_ylim, ylim_value + 2*h_ylim, ylim_value + h_ylim, ylim_value)) +
    theme(text=element_text(size=14)) +
    xlab("Method") + 
    ylab(ylables) + 
    ggtitle(paste0("Dropout Rate: ",dropout_rate,"%")) +
    scale_fill_discrete(name="Method",
                        breaks=c("1", "2", "3", "4", "5"),
                        labels=c("Dropout", "DrImpute", "scImpute", "MAGIC", "SCRABBLE")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  return(pp)
}

# calcualte the tsne
calculate_tsne <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # set up the parameter in the tsne calculation
  initial_dims_value = 30
  perplexity_value = 30
  
  # get the data
  data = get_data_HF(dropout_index, rand_seed)
  
  # defiine the list of tsne results
  set.seed(42)
  data_tsne = list()
  
  # calculate tsne of the true data
  data_tsne[[1]] = Rtsne(t(data$data_true), 
                         initial_dims = initial_dims_value, 
                         perplexity = perplexity_value)
  
  # calculate tsne of the raw data
  data_tsne[[2]] = Rtsne(t(data$data_raw), 
                         initial_dims = initial_dims_value, 
                         perplexity = perplexity_value)
  
  # calculate tsne of the imputed data using DrImpute
  data_tsne[[3]] = Rtsne(t(data$data_drimpute), 
                         initial_dims = initial_dims_value, 
                         perplexity = perplexity_value)
  
  # calculate tsne of the imputed data using scImpute
  data_tsne[[4]] = Rtsne(t(data$data_scimpute), 
                         initial_dims = initial_dims_value, 
                         perplexity = perplexity_value)
  
  # calculate tsne of the imputed data using MAGIC
  data_tsne[[5]] = Rtsne(t(data$data_magic), 
                         initial_dims = initial_dims_value, 
                         perplexity = perplexity_value)
  
  # calculate tsne of the imputed data using SCRABBLE
  data_tsne[[6]] = Rtsne(t(data$data_scrabble), 
                         initial_dims = initial_dims_value, 
                         perplexity = perplexity_value)
  
  # return the tsne result
  return(data_tsne)
  
}

# plot the tsne 
plot_tsne <- function(data, name, dot.size = 1){
  
  # Parameter in the function
  # data: the plot data consisting of three columns, the first 
  # two columns are the x and y of the tsne. The third one is 
  # the group information of the cells
  # name : the name of the plot
  # dot.size : the size of the dots
  
  # plot the figure
  p = ggplot(data, aes(x = x, y = y, color = ident)) +
    ggtitle(name) +
    xlab("TSNE_1") +
    ylab("TSNE_2") +
    geom_point(size = dot.size) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 12),
          legend.position="bottom",
          legend.title=element_blank()) +
    guides(colour=guide_legend(nrow=3,byrow=TRUE))
  
  return(p)
}

# plot tsne of the data
plot_tsne_HF <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # drop_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # calculate the tsne
  data_tsne = calculate_tsne(dropout_index, rand_seed)
  
  # define the list
  p = list()
  
  # define the names of the methods
  methods = c("True Data","Raw Data","DrImpute","scImpute","MAGIC","SCRABBLE")
  
  # plot the tsne figures
  for(i in c(1:6)){
    
    tmp = data_tsne[[i]]
    
    plot_data = as.matrix(tmp$Y)
    
    plot_data = data.frame(cbind(plot_data,rep(1:8,each = 100)))
    
    colnames(plot_data) = c("x","y","ident")
    
    p[[i]] = plot_tsne(plot_data, methods[i], dot.size = 1)
    
  }
  
  # gather the figures
  p1 = grid.arrange(grobs = p,ncol = 6)
  
  # return the figure handle
  return(p1)
}

# plot MA
plot_ma_HF <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # get the data
  data = get_data_HF(dropout_index, rand_seed)
  
  # calculat the dropout rate
  ratio_zeros_raw <- 1 - nnzero(data$data_raw)/length(data$data_raw)
  
  # define the list of figures
  pl = list()
  
  # plot MA of raw data
  pl[[1]] = prepare_ma_data_plot(data$data_true, data$data_raw, 
                                 paste0("Raw Data: ",round(100*ratio_zeros_raw),"%"))
  
  # plot MA of the imputed data of DrImpute
  pl[[2]] = prepare_ma_data_plot(data$data_true, data$data_drimpute, "DrImpute") 
  
  # plot MA of the imputed data of scImpute
  pl[[3]] = prepare_ma_data_plot(data$data_true, data$data_scimpute, "scImpute") 
  
  # plot MA of the imputed data of MAGIC
  pl[[4]] = prepare_ma_data_plot(data$data_true, data$data_magic, "MAGIC") 
  
  # plot MA of the imputed data of SCRABBLE
  pl[[5]] = prepare_ma_data_plot(data$data_true, data$data_scrabble, "SCRABBLE")
  
  # assemble of figues
  main = grid.arrange(grobs = pl,ncol = 5)
  
  # return the handle of the figure
  return(main)
}

prepare_ma_data_plot <- function(data_base,
                                 data,
                                 name_plot){
  
  # Paraemter in the function
  # data_base: the control data
  # data: the data
  # name_plot: the name of the plot
  
  # get the normalized data
  data = 10^data
  data_base = 10^data_base
  
  # construct the gene names
  genes = paste0("gene",c(1:dim(data_base)[1]))
  
  # construct the ma data.frame for the plots
  data_ma = data.frame(name = genes)
  
  data_ma$baseMean = rowMeans((log2(data_base) + log2(data))/2)
  
  data_ma$log2FoldChange = rowMeans(log2(data/data_base))
  
  data_ma$padj = t(calculate_pvalue_two_matrices(data_base,data))
  
  # MA plot
  p = ggmaplot1(data_ma,
                size = 1.5,
                top = 0,
                main = name_plot,
                font.main = c("bold",18),
                font.legend = c("bold", 15),
                legend = c(0.8, 0.2),
                genenames = as.vector(data_ma$name),
                ylim = c(-5,5)) +
    theme(axis.title = element_text(family = "Helvetica", size = (18)),
          axis.text.x = element_text(family = "Helvetica", size = (15)),
          axis.text.y = element_text(family = "Helvetica", size = (15)))
  
  
  return(p)
  
}



ggmaplot1 <- function (data, fdr = 0.05, fc = 1.5, genenames = NULL,
                       detection_call = NULL, size = NULL,
                       font.label = c(12, "plain", "black"), label.rectangle = FALSE,
                       palette = c("#B31B21", "#1465AC", "darkgray"),
                       top = 15, select.top.method = c("padj", "fc"),
                       main = NULL, xlab = "Log2 mean expression",  ylab = "Log2 fold change",
                       ggtheme = theme_classic(),...)
{
  
  if(!base::inherits(data, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
    stop("data must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if(!is.null(detection_call)){
    if(nrow(data)!=length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(data)")
  }
  else if("detection_call" %in% colnames(data)){
    detection_call = as.vector(data$detection_call)
  }
  else detection_call = rep(1, nrow(data))
  
  # Legend position
  if(is.null(list(...)$legend)) legend = c(0.12, 0.9)
  
  # Check data format
  ss = base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(data))
  if(length(ss)>0) stop("The colnames of data must contain: ",
                        paste(ss, collapse = ", "))
  
  if(is.null(genenames)) genenames <- rownames(data)
  else if(length(genenames)!=nrow(data))
    stop("genenames should be of length nrow(data).")
  
  sig = rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 1
  data = data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange,
                     padj = data$padj, sig = sig)
  
  # Change level labels
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- .levels(data$sig) %>% as.numeric()
  palette <- palette[.lev]
  new.levels <- c(
    paste0("Up: ", sum(sig == 1)),
    paste0("Down: ", sum(sig == 2)),
    "NS"
  ) %>% .[.lev]
  
  data$sig = factor(data$sig, labels = new.levels)
  
  
  # Ordering for selecting top gene
  select.top.method = match.arg(select.top.method)
  if(select.top.method == "padj") data <- data[order(data$padj), ]
  else if(select.top.method == "fc") data <- data[order(abs(data$lfc), decreasing = TRUE), ]
  # select data for top genes
  labs_data = stats::na.omit(data)
  labs_data = subset(labs_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))
  labs_data = utils::head(labs_data, top)
  
  font.label = .parse_font(font.label)
  font.label$size = ifelse(is.null(font.label$size), 12, font.label$size)
  font.label$color = ifelse(is.null(font.label$color), "black", font.label$color)
  font.label$face = ifelse(is.null(font.label$face), "plain", font.label$face)
  
  # Plot
  set.seed(42)
  mean <- lfc <- sig <- name <- padj <-  NULL
  p = ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig), size = size)
  
  if(label.rectangle){
    p = p + ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name),
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       force = 1, fontface = font.label$face,
                                       size = font.label$size/3, color = font.label$color)
  }
  else{
    p = p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name),
                                      box.padding = unit(0.35, "lines"),
                                      point.padding = unit(0.3, "lines"),
                                      force = 1, fontface = font.label$face,
                                      size = font.label$size/3, color = font.label$color)
  }
  
  p = p + scale_x_continuous(breaks=seq(0, max(data$mean), 2))+
    labs(x = xlab, y = ylab, title = main, color = "")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black"))
  
  p = ggpar(p, palette = palette, ggtheme = ggtheme, ...)
  
  # return the figure
  p
}

.levels <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  levels(x)
}

.parse_font <- function(font){
  if(is.null(font)) res <- NULL
  else if(inherits(font, "list")) res <- font
  else{
    # matching size and face
    size <- grep("^[0-9]+$", font, perl = TRUE)
    face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
    if(length(size) == 0) size <- NULL else size <- as.numeric(font[size])
    if(length(face) == 0) face <- NULL else face <- font[face]
    color <- setdiff(font, c(size, face))
    if(length(color) == 0) color <- NULL
    res <- list(size=size, face = face, color = color)
  }
  res
}

# plot the mean_variance
plot_meansd <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load all data
  data = get_data_HF(dropout_index, rand_seed)
  
  # calculate the dropout rate
  ratio_zeros_raw <- 1 - nnzero(data$data_raw)/length(data$data_raw)
  
  # define the list of figures
  pl = list()
  
  # plot the mean-variance of the true data
  p = meanSdPlot(data$data_true, ranks = FALSE)
  pl[[1]] = p$gg + ggtitle(paste0("True Data: ",round(100*ratio_zeros_true),"%")) 
  
  # plot the mean-variance of the raw data
  p = meanSdPlot(data$data_raw, ranks = FALSE)
  pl[[2]] = p$gg + ggtitle(paste0("Raw Data: ",round(100*ratio_zeros_raw),"%"))
  
  # plot the mean-variance of the imputed data of DrImpute
  p = meanSdPlot(data$data_drimpute, ranks = FALSE) 
  pl[[3]] = p$gg + ggtitle("DrImpute") 
  
  # plot the mean-variance of the imputed data of scImpute
  p = meanSdPlot(data$data_scimpute, ranks = FALSE)
  pl[[4]] = p$gg + ggtitle("scImpute")
  
  # plot the mean-variance of the imputed data of MAGIC
  p = meanSdPlot(data$data_magic, ranks = FALSE)
  pl[[5]] = p$gg + ggtitle("MAGIC")
  
  # plot the mean-variance of the imputed data of SCRABBLE
  p = meanSdPlot(data$data_scrabble, ranks = FALSE)
  pl[[6]] = p$gg + ggtitle("SCRABBLE")
  
  # Assemble the figures
  main = grid.arrange(grobs = pl,ncol = 6)
  
  # return the handler of figure
  return(main)
  
}

# plot the data
plot_data_p <- function(data, name){
  
  # Parameter in the function
  # data: the matrix, the column is the cell and the row is the gene
  # name: the name of the data
  
  
  # set the limitation of the data
  limit = c(0,5)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  colnames(data) = NULL
  rownames(data) = NULL
  
  # prepare the data for plotting
  longData = melt(as.matrix(data))
  pl = ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_colour_gradient2(limits=c(0, 5)) + 
    scale_fill_gradientn(colours = c("white", "blue", "red"), values = c(0,0.2,1)) +
    theme_bw()  + 
    scale_y_discrete(name ="Genes") +
    ggtitle(name) + 
    scale_x_discrete(name ="Cells") + 
    theme(panel.grid.major = element_blank(),
          legend.position="bottom",
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          line = element_blank(), 
          plot.title = element_text(family = "Helvetica", face = "bold", size = (8)),
          axis.title = element_text(family = "Helvetica", size = (6)),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    theme(legend.text=element_text(size=6),legend.title = element_text(size = 6))
  
  return(pl)
}

# plot the data
plot_data <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load all data
  data = get_data_HF(lambda_index, rand_seed, k)
  
  # calculate the dropout rate
  ratio_zeros_raw = 1 - nnzero(data$data_raw)/length(data$data_raw)
  
  # define the list of figures
  pl = list()
  
  # plot the true data
  pl[[1]] = plot_data_p(data$data_true, 
                        paste0("True Data: ",round(100*ratio_zeros_true),"%"))
  
  # plot the raw data
  pl[[2]] = plot_data_p(data$data_raw, 
                        paste0("Raw Data: ",round(100*ratio_zeros_raw),"%"))
  
  # plot the imputed data of DrImpute
  pl[[3]] = plot_data_p(data$data_drimpute, "DrImpute") 
  
  # plot the imputed data of scImpute
  pl[[4]] = plot_data_p(data$data_scimpute, "scImpute") 
  
  # plot the imputed data of MAGIC
  pl[[5]] = plot_data_p(data$data_magic, "MAGIC") 
  
  # plot the imputed data of SCRABBLE
  pl[[6]] = plot_data_p(data$data_scrabble, "SCRABBLE")
  
  # Assemble the figures
  main = grid.arrange(grobs = pl,ncol = 6)
  
  # return the handler of figure
  return(main)
  
}

# plot the correlation
plot_cor_p <- function(data, name){
  
  # data: the correlation matrix
  # name: the name of the data
  #
  
  # set up the data
  data1 = data
  diag(data1) = mean(c(data1))
  
  # setup the limitation
  limit = c(min(c(data1)),max(c(data1)))
  
  myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")))
  colnames(data) = NULL
  rownames(data) = NULL

  # prepare the plot data
  longData = melt(as.matrix(data))
  
  # plot the figure
  pl = ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    xlim(c(min(longData$Var1),max(longData$Var1))) +
    ylim(c(min(longData$Var2),max(longData$Var2))) +
    scale_colour_gradient2(limits=c(limit[1], limit[2])) + 
    scale_fill_gradientn(colours = c("blue", "white",  "red"), values = c(0,0.5,1)) +
    theme_bw()  + 
    scale_y_discrete(name ="Genes") +
    ggtitle(name) + 
    scale_x_discrete(name ="Cells") + 
    theme(panel.grid.major = element_blank(),
          legend.position="bottom",
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          line = element_blank(), 
          plot.title = element_text(family = "Helvetica", face = "bold", size = (8)),
          axis.title = element_text(family = "Helvetica", size = (6)),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    coord_fixed(ratio = 1) +
    theme(legend.text=element_text(size=6),legend.title = element_text(size = 6))
  
  # return the figure
  return(pl)
  
}

# plot the correlation data
plot_cor_HF <- function(dropout_index, rand_seed){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  # load all data
  data = get_data_HF(dropout_index, rand_seed)
  
  # set up the data
  data_true  = data$data_true
  data_dropout = data$data_raw
  data_drimpute = data$data_drimpute
  data_scimpute = data$data_scimpute
  data_magic = data$data_magic
  data_scrabble = data$data_scrabble
  
  # gene-gene correlation
  # true data
  data_true_gene = cor(t(data_true), method = "pearson")
  data_true_gene[is.na(data_true_gene)] = 0
  
  # dropout data
  data_dropout_gene = cor(t(data_dropout), method = "pearson")
  data_dropout_gene[is.na(data_dropout_gene)] = 0
  
  # DrImpute data
  data_drimpute_gene = cor(t(data_drimpute), method = "pearson")
  data_drimpute_gene[is.na(data_drimpute_gene)] = 0
  
  # scImpute data
  data_scimpute_gene = cor(t(data_scimpute), method = "pearson")
  data_scimpute_gene[is.na(data_scimpute_gene)] = 0
  
  # MAGIC data
  data_magic_gene = cor(t(data_magic), method = "pearson")
  data_magic_gene[is.na(data_magic_gene)] = 0
  
  # SCRABBLE data
  data_scrabble_gene = cor(t(data_scrabble), method = "pearson")
  data_scrabble_gene[is.na(data_scrabble_gene)] = 0
  
  p <- list()
  pl <- list()
  
  pl[[1]] <- plot_cor_p(data_true_gene,"Gene: True Data")
  
  pl[[2]] <- plot_cor_p(data_dropout_gene,"Gene: Dropout Data")
  
  pl[[3]] <- plot_cor_p(data_drimpute_gene,"Gene: DrImpute Data")
  
  pl[[4]] <- plot_cor_p(data_scimpute_gene,"Gene: scImpute Data")
  
  pl[[5]] <- plot_cor_p(data_magic_gene,"Gene: MAGIC Data")
  
  pl[[6]] <- plot_cor_p(data_scrabble_gene,"Gene: SCRABBLE Data")
  
  p[[1]] <- grid.arrange(grobs = pl,ncol = 6)
  
  # cell-cell correlation
  # true data
  data_true_cell = cor(data_true, method = "pearson")
  data_true_cell[is.na(data_true_cell)] = 0
  
  # dropout data
  data_dropout_cell = cor(data_dropout, method = "pearson")
  data_dropout_cell[is.na(data_dropout_cell)] = 0
  
  # DrImpute data
  data_drimpute_cell = cor(data_drimpute, method = "pearson")
  data_drimpute_cell[is.na(data_drimpute_cell)] = 0
  
  # scImpute data
  data_scimpute_cell = cor(data_scimpute, method = "pearson")
  data_scimpute_cell[is.na(data_scimpute_cell)] = 0
  
  # MAGIC data
  data_magic_cell = cor(data_magic, method = "pearson")
  data_magic_cell[is.na(data_magic_cell)] = 0
  
  # SCRABBLE
  data_scrabble_cell = cor(data_scrabble, method = "pearson")
  data_scrabble_cell[is.na(data_scrabble_cell)] = 0
  
  pl <- list()
  
  pl[[1]] <- plot_cor_p(data_true_cell,"Cell: True Data")
  
  pl[[2]] <- plot_cor_p(data_dropout_cell,"Cell: Dropout Data")
  
  pl[[3]] <- plot_cor_p(data_drimpute_cell,"Cell: DrImpute Data")
  
  pl[[4]] <- plot_cor_p(data_scimpute_cell,"Cell: scImpute Data")
  
  pl[[5]] <- plot_cor_p(data_magic_cell,"Cell: MAGIC Data")
  
  pl[[6]] <- plot_cor_p(data_scrabble_cell,"Cell: SCRABBLE Data")
  
  p[[2]] <- grid.arrange(grobs = pl,ncol = 6)
  
  return(p)
  
}














