# extract the data from bulk RNAseq data
get_filtered_mRNA_bulk <- function(tissue_name){
  
  # Parameter in the function
  # renName: the name of the tissue
  
  data_raw = read.table(file = paste0("data_raw/bulkRNAseq/19-tissues-expr/",tissue_name), 
                         header = T)
  
  # get the data
  data = data_raw[,c(1,6)]
  
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  mouse_id = getBM(attributes=c('mgi_symbol','ensembl_transcript_id','refseq_mrna'), 
                    filters = 'refseq_mrna', 
                    values = data$gene_id, 
                    mart = mouse)
  
  # get the gene name and gene ID
  data$gene_name = mouse_id$mgi_symbol[match(data$gene_id,mouse_id$refseq_mrna)]
  
  data$gene_ID = mouse_id$ensembl_transcript_id[match(data$gene_id,mouse_id$refseq_mrna)]
  
  # select the mRNA only
  data1 = data[!is.na(data$gene_name),]
  
  return(data1)
  
}

# get the data
get_data <- function(tissue_name,
                     data_tissue_file){
  
  # Parameter in the function
  # tissue_name: prepare the data of the tissue tissue_name
  # data_tissue_file: the file of the tissue in the bulk data
  
  # load all the scRNAseq data
  load(file = paste0("data_sc_bulk/sc_",tissue_name,".RData"))
  
  # get the file name of the bulk RNAseq data 
  bulk_file = data_tissue_file$V1[data_tissue_file$V2 == tissue_name]
  
  # define the list for the bulk RNAseq data
  bulk_data_tmp = list()
  
  # the gene list of the scRNAseq data
  common_gene = rownames(data_select)
  
  # get the bulk RNAseq data
  for(i in c(1:length(bulk_file))){
    
    tmp = get_filtered_mRNA_bulk(bulk_file[i])
    
    common_gene = intersect(common_gene,tmp$gene_name)
    
    bulk_data_tmp[[i]] = tmp
    
  }
  
  # combine bulk RNAseq data
  data_bulk = c()
  
  # make the gene list of each data consistent
  for(i in c(1:length(bulk_file))){
  
      tmp = bulk_data_tmp[[i]] 
    
      data_bulk = cbind(data_bulk, tmp$FPKM[match(common_gene,tmp$gene_name)])
      
  }
  
  # get the bulk RNAseq for the imputation
  data_bulk_avg = log1p(rowMeans(data_bulk))
  
  data_bulk_avg = as.data.frame(data_bulk_avg)
  
  rownames(data_bulk_avg) = common_gene
  
  # get the scRNAseq data
  index_sc = match(common_gene, rownames(data_select))
  
  data_select1 = as.data.frame(data_select)
  
  data_select1 = data_select1[index_sc,]
  
  rownames(data_select1) = common_gene
  
  # prepare the Seurat object
  tiss = CreateSeuratObject(raw.data = data_select1)
  
  tiss = AddMetaData(object = tiss, meta_data)
  
  tiss = NormalizeData(object = tiss)
  
  tiss = ScaleData(object = tiss)
  
  tiss = FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.0001)
  
  data1 = GetAssayData(tiss,slot = "data")
  
  data1 = as(data1,"matrix")
  
  # normalize the scRNAseq and bulk RNAseq data
  tmp = MatNonZeroMeanSd(data1)
  
  ratio = data_bulk_avg$data_bulk_avg[which.max(tmp$nzeros)]/mean(data1[which.max(tmp$nzeros),])
  
  data_bulk_avg$data_bulk_avg = data_bulk_avg$data_bulk_avg/ratio
  
  index_sc = match(tiss@var.genes,common_gene)
  
  # define the imputation list
  data = list()
  
  data[[1]] = data1[index_sc,]
  
  data[[2]] = data_bulk_avg[index_sc,]
  
  data[[3]] = tiss@var.genes
  
  data[[4]] = data1
  
  data[[5]] = data_bulk_avg
  
  # save the data
  saveRDS(data,file = paste0("data_sc_bulk/",tissue_name,"_imputation.rds"))
  
}

# plot tsne
plot_tsne <- function(plot_data,
                      name, 
                      dot.size = 1){
  # Parameter in the function
  # plot_data: the data for visualization
  # name: the name of the plot
  # do.size: the size of the dot
  
  plot_data %>%
    dplyr::group_by(ident) %>%
    summarize(x = median(x = x), y = median(x = y)) -> centers
  
  # plot the dots
  p = ggplot(plot_data,aes(x = x, y = y, color = as.factor(ident))) +
    ggtitle(name) +
    xlab("TSNE_1") +
    ylab("TSNE_2") +
    geom_point(size = dot.size) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 18, hjust = 0.4),
          axis.text = element_text(size = 12),
          legend.title=element_blank())
  
  # plot the annotation number
  p = p + geom_text(data = centers, mapping = aes(x = x, y = y, label = ident), colour = "black")
  
  
  return(p)
  
}


calculate_cluster <- function(cluster_index,
                             cluster,
                             data_tsne){
  
  # Parameter in the function
  # cluster_index: the index of the cluster used for calculating dun index
  # cluster: the clustering informaiton of the cells
  # data_tsne: the tsne data

  cluster1 = cluster
  
  index = cluster1 != cluster_index
  
  cluster1[index] = max(cluster) + 1
  
  stats_cluster_dun = matrix(0,nrow = 5,ncol = 1)
  
  # calculate the tsne
  stat_tmp = cluster.stats(dist(as.matrix(data_tsne[[1]]$Y)),cluster1)
  
  stats_cluster_dun[1] = stat_tmp$dunn

  stat_tmp = cluster.stats(dist(as.matrix(data_tsne[[2]]$Y)),cluster1)
  
  stats_cluster_dun[2] = stat_tmp$dunn

  stat_tmp = cluster.stats(dist(as.matrix(data_tsne[[3]]$Y)),cluster1)
  
  stats_cluster_dun[3] = stat_tmp$dunn

  stat_tmp = cluster.stats(dist(as.matrix(data_tsne[[4]]$Y)),cluster1)
  
  stats_cluster_dun[4] = stat_tmp$dunn

  stat_tmp = cluster.stats(dist(as.matrix(data_tsne[[5]]$Y)),cluster1)
  
  stats_cluster_dun[5] = stat_tmp$dunn

  return(stats_cluster)
  
}

# plot the bun values
plot_bun_value <- function(bun_value,annotation_info){
  
  # Parameter in the function
  # bun_value: the bun values
  # annotation_info: the cell annotation information
  
  # calculate the log2 values
  index_tmp = -log2(bun_value)
  
  # prepare the data for the boxplot
  longData = melt(as.matrix(index_tmp))
  
  # define the column name
  colnames(longData) = c("X1","X2","value")
  
  # plot the boxplot
  pp = ggplot(data=longData, aes(x=as.factor(X2), y= value, fill= as.factor(X1))) +
    geom_bar(stat="identity", position=position_dodge()) + 
    scale_color_manual(labels = c("Raw", "DrImpute","scImpute","MAGIC","SCRABBLE")) +
    scale_fill_manual(values= c("#00AFBB","#0000CD","#E7B800", "#FC4E07", "#6ebb00"), 
                      name="Data Type",
                      breaks=c(1,2,3,4,5),
                      labels=c("Raw Data","DrImpute","scImpute","MAGIC","SCRABBLE")) +
    theme_bw() +
    ylab("-log2(Dunn Index)") +
    scale_x_discrete(breaks= annotation_info[,2],labels=annotation_info[,1]) +
    theme(axis.title.x=element_blank(),axis.text.x = element_text(size=9, angle=45 ,vjust=0.6),
          axis.title.y=element_text(size=14),axis.text.y = element_text(size=12),
          legend.position = "bottom") +
    guides(fill = guide_legend(title = "Data Type"))
  
  return(pp)
  
}

pdf_dun_tsne <- function(tissue_name){
  
  # Parameter in the function
  # tissue_name: the name of the tissue
  
  # load the tSNE data
  data_tsne = readRDS(file = paste0("data_all/",tissue_name,"_TSNE.rds"))
  
  # load the scRNAseq and related data
  load(paste0("data_sc_bulk/sc_",tissue_name,".RData"))
  
  pl = list()
  
  meta_data$Cluster = t(data.frame(strsplit(meta_data$ClusterID, "_"))[2,])[,1]
  
  tmp = as.numeric(meta_data$Cluster)
  
  tmp[is.na(tmp)] = 0
  
  knn5 = knn(data_tsne[[5]]$Y, data_tsne[[5]]$Y, tmp, k=5)
  
  knn5 = as.numeric(knn5)
  
  method_name = c("Raw Data","DrImpute","scImpute","MAGIC","SCRABBLE")
  
  for(i in c(1:5)){
    
    tmp = data_tsne[[i]]
    
    plot_data = as.matrix(tmp$Y)
    
    plot_data = data.frame(cbind(plot_data,as.numeric(knn5)))
    
    colnames(plot_data) = c("x","y","ident")
    
    pl[[i]] = plot_tsne(plot_data,paste0(method_name[i],": ",tissue_name), dot.size = 1)
    
  }
  
  p1 = grid.arrange(grobs = pl,ncol = 3)
  
  ggsave(
    filename = paste0(tissue_name, "_tsne.pdf"),
    plot = p1,
    width = 20,
    height = 10
  )
  
  annotation_cell = cbind(meta_data$Annotation,meta_data$Cluster)
  
  index = duplicated(as.data.frame(annotation_cell[,2]))
  
  annotation_info = annotation_cell[!index,]
  
  annotation_info = annotation_info[order(as.numeric(annotation_info[,2])),]
  
  cluster_n = unique(knn5)
  
  bun_value = c()
  
  for(i in c(1:dim(annotation_info)[1])){
    
    if(i %in% cluster_n){
      
      tmp_val = calculate_cluster(i,as.numeric(knn5),data_tsne)
      
      bun_value = cbind(bun_value,tmp_val[[1]])
      
    }else{
      
      bun_value = cbind(bun_value,c(1,1,1,1,1))
      
    }
    
  }
  
  pp = plot_bun_value(bun_value,annotation_info)
  
  ggsave(
    filename = paste0(tissue_name, "_tsne_bun_index.pdf"),
    plot = pp,
    width = 8,
    height = 6
  )
}


pdf_dun_tsne1 = function(tissue_name){
  
  # Parameter in the function
  # tissue_name: the name of the tissue
  
  # load the tSNE data
  data_tsne = readRDS(file = paste0("data_all/",tissue_name,"_TSNE.rds"))
  
  # load the scRNAseq and related data
  load(paste0("data_sc_bulk/sc_",tissue_name,".RData"))
  
  pl = list()
  
  meta_data$Cluster = t(data.frame(strsplit(meta_data$ClusterID, "_"))[3,])[,1]
  
  tmp = as.numeric(meta_data$Cluster)
  
  tmp[is.na(tmp)] = 0
  
  knn5 = knn(data_tsne[[5]]$Y, data_tsne[[5]]$Y, tmp, k=5)
  
  method_name = c("Raw Data","DrImpute","scImpute","MAGIC","SCRABBLE")
  
  for(i in c(1:5)){
    
    tmp = data_tsne[[i]]
    
    plot_data = as.matrix(tmp$Y)
    
    plot_data = data.frame(cbind(plot_data,as.numeric(knn5)))
    
    colnames(plot_data) = c("x","y","ident")
    
    pl[[i]] = plot_tsne(plot_data,paste0(method_name[i],": ",tissue_name), dot.size = 1)
    
  }
  
  p1 = grid.arrange(grobs = pl,ncol = 3)
  
  ggsave(
    filename = paste0(tissue_name, "_tsne.pdf"),
    plot = p1,
    width = 20,
    height = 10
    
  )
  
  annotation_cell = cbind(meta_data$Annotation,meta_data$Cluster)
  
  index = duplicated(as.data.frame(annotation_cell[,2]))
  
  annotation_info = annotation_cell[!index,]
  
  annotation_info = annotation_info[order(as.numeric(annotation_info[,2])),]
  
  cluster_n = unique(knn5)
  
  bun_value = c()
  
  for(i in c(1:dim(annotation_info)[1])){
    
    if(i %in% cluster_n){
      
      tmp_val = calculate_cluster(i,as.numeric(knn5),data_tsne)
      
      bun_value = cbind(bun_value,tmp_val[[1]])
      
    }else{
      
      bun_value = cbind(bun_value,c(1,1,1,1,1))
      
    }
    
  }
  
  pp = plot_bun_value(bun_value,annotation_info)
  
  ggsave(
    filename = paste0(tissue_name, "_tsne_bun_index.pdf"),
    plot = pp,
    width = 8,
    height = 6
  )
  
}



