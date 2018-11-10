# plot the figures of data
plot_data <- function(data,name){
  limit <- c(0,4)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  colnames(data) <- NULL
  rownames(data) <- NULL
  longData<-melt(as.matrix(data))
  pl <- ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_colour_gradient2(limits=c(0, 4)) +
    scale_fill_gradient2(limits=c(0, 4),low = "blue", mid = "white",high = "red", midpoint = 2) +
    theme_bw()  +
    scale_y_discrete(name ="Genes") +
    ggtitle(name) +
    scale_x_discrete(name ="Cells") +
    theme(panel.grid.major = element_blank(),
          legend.position="bottom",
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          line = element_blank(),
          plot.title = element_text(family = "Helvetica", face = "bold", size = (12)),
          axis.title = element_text(family = "Helvetica", size = (10)),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    theme(legend.text=element_text(size=10),legend.title = element_text(size = 10))
  return(pl)
}
