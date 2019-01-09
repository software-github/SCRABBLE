# repeat the rows to make a matrix
rep.row <- function(x,n){
  
  matrix(rep(x,each=n),nrow=n)
  
}

# repeat the columns to make a matrix
rep.col <- function(x,n){
  
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  
}

# plot the distributions
distribution_FISH <- function(data_sc, 
                              data_drimpute, 
                              data_scimpute, 
                              data_magic, 
                              data_scrabble, 
                              gene_name, 
                              gene1,
                              ratio) {
  
  # Parameter in the function
  # data_sc: the single cell data
  # data_drimpute: the imputed data of drimpute
  # data_scimpute: the imputed data of scimpute
  # data_magic: the imputed data of magic
  # data_scrabble: the imputed data of scrabble
  # gene_name: the gene list of the scRNAseq data
  # gene1: the gene whose distribution plot is drawn
  # ratio: the normalized constant between the scRNAseq and smFISH data
  
  # get the index of the gene
  index = which(gene_name$V1 == gene1)
  
  # get raw data
  tmp = data.frame(t(data_sc[index,]))
  
  colnames(tmp) = "value"
  
  tmp1 = log10(tmp+1)
  
  tmp1$e = 1
  
  # get drimpute data
  tmp = data.frame((data_drimpute[index,]))
  
  colnames(tmp) = "value"
  
  tmp2 = log10(tmp+1)
  
  tmp2$e = 2
  
  # get scimpute data
  tmp = data.frame(t(data_scimpute[index,]))
  
  colnames(tmp) = "value"
  
  tmp3 = log10(tmp+1)
  
  tmp3$e = 3
  
  # get magic data
  tmp = data.frame(t(data_magic[index,])) 
  
  colnames(tmp) = "value"
  
  tmp4 = log10(tmp+1)
  
  tmp4$e = 4
  
  # get scrabble data
  tmp = data.frame((data_scrabble[index,]))
  
  colnames(tmp) = "value"
  
  tmp5 = log10(tmp+1)
  
  tmp5$e = 5
  
  # get smFISH data
  data_FISH = read.table(file = paste0("data_mFISH/",gene1,".csv"), 
                         header = F, 
                         sep=",", 
                         stringsAsFactors = F)
  
  tmp6 = ratio*log10(data_FISH+1)
  
  colnames(tmp6) = "value"

  tmp6$e <- 6
  
  # combine the data
  data = rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
  
  colnames(data) = c("value","e")
  
  # plot the density graph
  p = ggplot(data, aes(x = value,color = as.factor(e), fill=as.factor(e))) +
    geom_density( alpha = 0.2,adjust = 15) +
    theme_bw() + 
    scale_fill_manual(labels = c("Dropseq","DrImpute","scImpute","MAGIC","SCRABBLE","smFISH"),
                      values=c("#00AFBB","#0000CD", "#E7B800", "#FC4E07", "#6ebb00","#FF1493")) + 
    scale_color_manual(labels = c("Dropseq","DrImpute", "scImpute","MAGIC","SCRABBLE","smFISH"),
                       values=c("#00AFBB","#0000CD", "#E7B800","#FC4E07",  "#6ebb00","#FF1493")) +
    ggtitle(paste0(gene1)) +
    theme(legend.position = "bottom",
          legend.title=element_blank()) +
    labs(x="log10(Normalized Expression levels + 1)", y="Density")+
    scale_x_continuous(limits=c(-1,5)) +
    scale_y_continuous(limits = c(0,1))
  
  # return the handler of plot
  return(p)
  
}