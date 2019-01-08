performance_comparison <- function(gene_only_SCRB,
                                   gene_filter,
                                   threshold,
                                   data_de){
  
  # Parameter in the function
  # gene_only_SCRB: the gene list without dropouts in SCRB
  # gene_filter: the gene list of the data
  # threshold: the threshold of genes
  # data_de: the combined data
  
  
  # get the DropSeq RNAseq Data
  data_sc_Dropseq = data_de[,dataType == 2]
  
  # get genes whose dropouts are greater than threshold in DropSeq
  tmp_dropseq = data_sc_Dropseq[match(gene_only_SCRB,gene_filter),]
  
  index = rowSums(tmp_dropseq == 0) > threshold
  
  colnames(index) = NULL
  
  common_diff = gene_only_SCRB[which(index)]
  
  # get the data
  data_ks = data_de[match(common_diff,gene_filter),]
  
  # get SCRB data
  data_ks_SCRB = data_ks[,dataType == 1]
  
  # get DropSeq data
  data_ks_drop = data_ks[,dataType == 2]
  
  # gete the imputed data of DrImpute
  data_ks_drimpute = data_ks[,dataType == 3]
  
  # get the imputed data of scImpute
  data_ks_scimpute = data_ks[,dataType == 4]
  
  # get the imputed data of MAGIC
  data_ks_magic = data_ks[,dataType == 5]
  
  # get the imputed data of SCRABBLE
  data_ks_scrabble = data_ks[,dataType == 6]
  
  # define the KS statistics matrix
  ks_value = matrix(0,nrow = dim(data_ks)[1], ncol = 5)
  
  # Calculate the KS statistics
  for (i in c(1:dim(data_ks)[1])){
    
    # raw data vs SCRB
    tmp = ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_drop[i,]))
    
    ks_value[i,1] = tmp$statistic
    
    # DrImpute vs SCRB
    tmp = ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_drimpute[i,]))
    
    ks_value[i,2] = tmp$statistic
    
    # scImpute vs SCRB
    tmp = ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_scimpute[i,]))
    
    ks_value[i,3] = tmp$statistic
    
    # MAGIC vs SCRB
    tmp = ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_magic[i,]))
    
    ks_value[i,4] = tmp$statistic
    
    # SCRABBLE vs SCRB
    tmp = ks.test(as.matrix(data_ks_SCRB[i,]), as.matrix(data_ks_scrabble[i,]))
    
    ks_value[i,5] = tmp$statistic
    
  }
  
  # prepare the data for boxplot
  tmp = melt(ks_value)
  
  # set up as the data.frame
  ksData = data.frame(ks=tmp$value,group=tmp$Var2)
  
  # set up the comparison
  my_comparisons = list( c("1", "5"),c("2", "5"), c("3", "5"), c("4", "5"))
  
  # compute the pvalue
  pval1 = compare_means(ks ~ group,data = ksData, method = "t.test", paired = T)
  
  pval = c(paste0('p1 = ',formatC(pval1$p[4], format = "e", digits = 2)),
           paste0('p1 = ',formatC(pval1$p[7], format = "e", digits = 2)),
           paste0('p1 = ',formatC(pval1$p[9], format = "e", digits = 2)),
           paste0('p1 = ',formatC(pval1$p[10], format = "e", digits = 2)))
  
  # plot the boxplot
  pl = ggboxplot(ksData, x = "group", y = "ks", 
                 fill = "group",
                 palette = c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00")) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    ylim(c(0,1.5)) + 
    theme_bw() +
    ggtitle(paste0("Genes used: ", length(common_diff), 
                   "  Dropout cutoff: ", 
                   formatC(threshold*100/76, digits = 0, format = "f"),"%")) +
    geom_signif(comparisons = my_comparisons, 
                annotations = pval,
                tip_length = 0.03,
                y_position = c(1.4,1.3, 1.2, 1.1)) +
    theme(text=element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Data Sets") + 
    ylab("K-S Statistics")
  
  # return the handler of boxplot
  return(pl)
  
}

performance_distribution_comparison <- function(gene_only_SCRB,
                                                gene_filter,
                                                threshold,
                                                data_de){
  
  # Parameter in the function
  # gene_only_SCRB: the gene list without dropouts in SCRB
  # gene_filter: the gene list of the data
  # threshold: the threshold of genes
  # data_de: the combined data
  
  # get the DropSeq RNAseq Data
  data_sc_Dropseq = data_de[,dataType == 2]
  
  # get genes whose dropouts are greater than threshold in DropSeq
  index = rowSums(tmp_dropseq == 0) > threshold
  
  colnames(index) = NULL
  
  common_diff = gene_only_SCRB[which(index)]
  
  # get the data
  data_ks = data_de[match(common_diff,gene_filter),]
  
  # get SCRB data
  data_ks_SCRB = data_ks[,dataType == 1]
  
  # get DropSeq data
  data_ks_drop = data_ks[,dataType == 2]
  
  # gete the imputed data of DrImpute
  data_ks_drimpute = data_ks[,dataType == 3]
  
  # get the imputed data of scImpute
  data_ks_scimpute = data_ks[,dataType == 4]
  
  # get the imputed data of MAGIC
  data_ks_magic = data_ks[,dataType == 5]
  
  # get the imputed data of SCRABBLE
  data_ks_scrabble = data_ks[,dataType == 6]
  
  # define the list
  p = list()
  
  for (i in c(1:dim(data_ks)[1])){
    
    # extract the data of scrb
    tmp1 = data.frame( value = t(data_ks_SCRB[i,]))
    
    tmp1$e = 6 
    
    # extract the data of raw data
    tmp2 = data.frame( value = t(data_ks_drop[i,]))
    
    tmp2$e = 1
    
    # extract the data of drimpute imputed data
    tmp3 = data.frame( value = t(data_ks_drimpute[i,]))
    
    tmp3$e = 2
    
    # extract the data of scimpute imputed data
    tmp4 = data.frame( value = t(data_ks_scimpute[i,]))
    
    tmp4$e = 3
    
    # extract the data of magic imputed data
    tmp5 = data.frame( value = t(data_ks_magic[i,]))
    
    tmp5$e = 4
    
    # extract the data of scrabble imputed data
    tmp6 = data.frame( value = t(data_ks_scrabble[i,]))
    
    tmp6$e = 5
    
    # assemble all data as a dataframe
    data = rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
    
    colnames(data) = c("value","e")

    # plot the distribution density 
    p[[i]] = ggplot(data, aes(x = value,color = as.factor(e), fill=as.factor(e))) +
      geom_density( alpha = 0.2,adjust = 10) +
      theme_bw() + 
      scale_fill_manual(labels = c("Dropseq","DrImpute","scImpute","MAGIC","SCRABBLE","SCRBseq"),
                        values=c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00","#8B008B")) + 
      scale_color_manual(labels = c("Dropseq","DrImpute", "scImpute","MAGIC","SCRABBLE","SCRBseq"),
                         values=c("#00AFBB","#0000CD", "#E7B800","#FC4E07",  "#6ebb00","#8B008B")) +
      ggtitle(paste0(common_diff[i])) +
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      labs(x="log10(Normalized Expression levels + 1)", y="Density")+
      scale_x_continuous(limits=c(-1,5)) +
      scale_y_continuous(limits = c(0,1))
  }
  
  # return the handler of the plot
  return(p)
  
}
