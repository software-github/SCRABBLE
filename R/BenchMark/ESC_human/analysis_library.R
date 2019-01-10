# calculate the variance across the rows or the columns
MatVar <- function(x, dim = 1, ...) {
  
  # parameter in the function
  # x: the input matrix 
  # dim: determine the variance across the row or column.
  # 1 is the one across the row and 2 is the one across the column
  
  if(dim == 1){
    
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    
  } else if (dim == 2) {
    
    rowSums((t(x) - colMeans(x, ...))^2, ...)/(dim(x)[1] - 1)
    
  } else stop("Please enter valid dimension")
}

# plot the pathway ratio
plot_pathway_ratio <- function(dataV3, dataType, pathway_name){
  
  # Parameter in the function
  # dataV3: the ratio of the different methods, the columns are 
  # the different methods
  # dataType: the cell types
  # pathway_name: the name of the pathway database
  
  if(dataType == "NPC"){
    
    # define the number of data points
    N = dim(dataV3)[2]
    
    # define the lowest levels of the pvalues
    hlim = 190
    
    # prepare the data as a vector for the boxplot
    dataV1 = data.frame(y = (as.vector(as.matrix(dataV3))))
    
    dataV1$group = rep(c(1:5),N)
    
    # calculate the pvalues
    my_comparisons = list( c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))
    
    pval = compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "5")
    
    # plot the boxplots
    pl = ggboxplot(dataV1, x = "group", y = "y", fill = "group",
                    palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00"),outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", width = 0.3) + 
      ylim(c(-50,1.45*hlim)) + 
      theme_bw() +
      ggtitle(paste0(dataType,"_",pathway_name)) +
      geom_signif(comparisons = my_comparisons, 
                  annotations = formatC(pval$p, format = "e", digits = 2),
                  tip_length = 0.03,
                  y_position = c(1.4*hlim,1.3*hlim, 1.2*hlim, 1.1*hlim, hlim)) +
      theme(text=element_text(size=12),legend.position="bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      xlab("Data Sets") + 
      ylab("Increasing Percentage")
    
  }else{
    
    # define the number of data points
    N = dim(dataV3)[2]
    
    # define the lowest levels of the pvalues
    hlim = 150
    
    # prepare the data as a vector for the boxplot
    dataV1 = data.frame(y = (as.vector(as.matrix(dataV3))))
    
    # calculate the pvalues
    dataV1$group = rep(c(1:5),N)
    
    my_comparisons = list( c("1", "5"), c("2", "5"), c("3", "5"), c("4", "5"))
    
    pval = compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "5")
    
    # plot the boxplots
    pl = ggboxplot(dataV1, x = "group", y = "y", fill = "group",
                    palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00"),outlier.shape = NA) +
      stat_boxplot(geom = "errorbar", width = 0.3) + 
      ylim(c(-50,1.45*hlim)) + 
      theme_bw() +
      ggtitle(paste0(dataType,"_",pathway_name)) +
      geom_signif(comparisons = my_comparisons, 
                  annotations = formatC(pval$p, format = "e", digits = 2),
                  tip_length = 0.03,
                  y_position = c(1.4*hlim,1.3*hlim, 1.2*hlim, 1.1*hlim, hlim)) +
      theme(text=element_text(size=12),legend.position="bottom",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      xlab("Data Sets") + 
      ylab("Increasing Percentage")
    
  }
  
  return(pl)
  
}


generate_index <- function(data_set, pathway_name){
  
  # Parameter in the function
  # data_set: the data of cell types
  # pathway_name: the name of pathway database
  
  # load the gene
  data_gene = fread(file = "data_all/gene_name.csv", 
                    sep=',', header = FALSE)
  
  # load the scRNAseq data
  data_sc = as.matrix(fread(file = paste0("data_all/data_sc_",
                                           data_set,
                                           ".csv")))
  
  # calculate the variance across the genes
  var0 = MatVar(data_sc,1)
  
  index_sc = var0 > 1e-10
  
  data_gene = data_gene[index_sc,]
  
  n_gene = dim(data_gene)[1]
  
  # determine the pathway numbers
  if(pathway_name == "IPA"){
    
    N = 186
    
  }
  
  if(pathway_name == "KEGG"){
    
    N = 186
    
  }
  
  if(pathway_name == "REACTOME"){
    
    N = 674
    
  }
  
  # calculate the index
  index = list()
  
  for(i in c(1:N)){
    
    tmp = fread(file = paste0("data_all/", pathway_name,"/",pathway_name,"_gene_",i,".csv"), 
                 header = FALSE)
    
    tmp = match(tmp$V1, data_gene$V1)
    
    tmp = tmp[!is.na(tmp)]
    
    index[[i]] = tmp
  }
  
  for(i in c(1:N)){
    
    tmp = c()
    
    k = length(index[[i]])
    
    for(j in c(1:100)){
      
      set.seed(j)
      
      tmp = rbind(tmp, sample(1:n_gene,k))
    }
    
    index[[i]] = rbind(index[[i]],tmp)
    
  }
  
  save(index, index_sc, file = paste0("data_all/",
                                      pathway_name,
                                      "_",
                                      data_set,
                                      "_index.RData")
  )
} 

# Calculate the ratio
calculate_ratio <- function(data_set, pathway_name, method_name){
  
  # Parameter in the function
  # data_set: the data of cell type
  # pathway_name: the name of the pathway
  # method_name: the name of the method used for the imputation
  
  load(file = paste0("data_all/",
                     pathway_name,
                     "_",
                     data_set,
                     "_index.RData")
  )
  
  
  if(method_name == "dropout"){
    
    data_sc = as.matrix(fread(file = paste0("data_all/data_sc_",
                                             data_set,
                                             ".csv")))
  }
  
  if(method_name == "drimpute"){
    
    data_sc = readRDS(file = paste0("imputation_drimpute_data/data_",
                                     data_set,
                                     "_drimpute_imputation.rds"))
  }
  
  if(method_name == "scimpute"){
    
    data_sc = readRDS(file = paste0("imputation_scimpute_data/data_",
                                     data_set,
                                     "_scimpute_imputation.rds"))
  }
  
  if(method_name == "magic"){
    
    data_tmp = as.matrix(fread(file = paste0("imputation_magic_data/data_magic_",
                                              data_set,
                                              ".csv")))
    
    data_sc = data_tmp[,-1] 
  }
  
  if(method_name == "scrabble"){
    
    data_sc = readRDS(file = paste0("imputation_scrabble_data/data_",
                                     data_set,
                                     "_scrabble_imputation.rds"))
  }
  
  # determine the pathway numbers
  if(pathway_name == "IPA"){
    
    N = 186
    
  }
  
  if(pathway_name == "KEGG"){
    
    N = 186
    
  }
  
  if(pathway_name == "REACTOME"){
    
    N = 674
    
  }
  
  values = c()
  
  data_sc = data_sc[index_sc,]
  
  for(i in c(1:N)){
    
    tmp_index = index[[i]]
    
    N_index = dim(tmp_index)[2]
    
    mean_value = c()
    
    if (N_index > 10){
      
      for(j in c(1:101)){
        
        tmp_data = data_sc[tmp_index[j,],]
        
        cor_tmp = cor(t(tmp_data))
        
        tmp1 = abs(cor_tmp[lower.tri(cor_tmp)])
        
        tmp1 = tmp1[!is.nan(tmp1)]
        
        mean_value = cbind(mean_value, mean(tmp1))
        
      }
      
      values = rbind(values, c(mean_value[1], mean(mean_value[2:101])))
    }
  }
  
  values = values[!is.na(values[,2]),]
  
  return(values)
  
}
