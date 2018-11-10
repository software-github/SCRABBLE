# SCRABBLE

### SCRABBLE has been implemented in R and MATLAB.

SCRABBLE imputes drop-out data by optimizing an objective function that consists of three terms. The first term ensures that imputed values for genes with nonzero expression remain as close to their original values as possible, thus minimizing unwanted bias towards expressed genes. The second term ensures the rank of the imputed data matrix to be as small as possible. The rationale is that we only expect a limited number of distinct cell types in the samples. The third term operates on the bulk RNA-Seq data. It ensures consistency between the average gene expression of the aggregated imputed data and the average gene expression of the bulk RNA-Seq data. We developed a convex optimization algorithm to minimize the objective function.

### Quick tutorial
```
data_sc <- data[[1]]
data_bulk <- data[[2]]
data_true <- data[[3]]

parameter <- c(10,1e-5,1e-4)
nIter <- 20

result <- scrabble(data1,
                   parameter = parameter, 
                   nIter = 30,
                   error_out_threshold = 1e-7, 
                   nIter_inner = 100,
                   error_inner_threshold = 1e-5)
```                   
