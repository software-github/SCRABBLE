# plot the figures of data
#'
#' It is used to plot the data as a images
#'
#' @param data the input data which is a data matrix.
#'
#' @param name the name of the data and it is a string.
#'
#'
#' @examples
#' # Run plot_data
#' p <- plot_data(demo_data[[1]],"Drop-out Data")
#'
#' @return ggplot handler
#'
#' @export
#'
#'

plot_data <- function(data,name){
  limit <- c(0,5)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  print(dim(data))
  colnames(data) <- NULL
  rownames(data) <- NULL
  longData<-melt(as.matrix(data))
  colnames(longData) <- c("Var1", "Var2","value")
  pl <- ggplot(longData, aes_string(x = "Var2", y = "Var1")) +
    geom_raster(aes_string(fill= "value")) +
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
