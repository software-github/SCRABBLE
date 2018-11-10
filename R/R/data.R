#' Test data for scrabble
#'
#' @description  "data" is a data list with the length of 3. The first element in
#' the list is generated drop-out scRNAseq data with 732 genes and 1000 cells. The second
#' element in the list is the generated bulk RNAseq data with 732 genes. The third
#' element is the true scRNAseq data without dropouts. The steps of generating the data
#' is shown in Details section.
#'
#' @usage data_sc <- data[[1]]
#' data_bulk <- data[[2]]
#' data_true <- data[[3]]
#'
#' @author Tao Peng, Kai Tan
#'
#' @details The data set was generated from down sampling from bulk RNAseq data.
#' We used the bulk RNA-Seq data set of mouse hair follicles (GSE85039).
#' In total, the dataset contains 20 different combinations of anatomic
#' sites and developmental time points, thus constituting a high dimensional
#' measurement space. We used the following procedures to generate the
#' drop-out datasets. 1) We selected 732 genes that are
#' differentially expressed in the 20 conditions based on ANOVA analysis.
#' 2) We randomly selected 10 out of the 20 conditions.  3) For each condition,
#' we generated 100 resampled datasets. The means and standard deviations of
#' genes were calculated for each condition based on the 100 resampled datasets.
#' 4) 100 new datasets were generated based on the mean and the standard deviation
#' of each gene. 5) The final data set was obtained by combining 1000 samples
#' representing the 10 conditions. This 1000x732 matrix now represents 1000 cells
#' and 732 genes. 6) we make the drop-out rate of each gene in each cell following
#' a double exponential function . Zero values are introduced into the simulated
#' data for each gene in each cell based on the Bernoulli distribution defined by
#' the corresponding drop-out rate.
#'
"data"
