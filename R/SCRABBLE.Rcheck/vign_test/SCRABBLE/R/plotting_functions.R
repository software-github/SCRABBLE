# plot the figures of data
plot_data <- function(data,name){
  limit <- c(0,5)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  print(dim(data))
  colnames(data) <- NULL
  rownames(data) <- NULL
  longData<-melt(as.matrix(data))
  pl <- ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_colour_gradient2(limits=c(0, 5)) +
    scale_fill_gradientn(colours = c("white", "blue", "red"), values = c(0,0.6,1)) +
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
