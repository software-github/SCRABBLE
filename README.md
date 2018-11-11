# SCRABBLE

## SCRABBLE has been implemented in R and MATLAB.

SCRABBLE imputes drop-out data by optimizing an objective function that consists of three terms. The first term ensures that imputed values for genes with nonzero expression remain as close to their original values as possible, thus minimizing unwanted bias towards expressed genes. The second term ensures the rank of the imputed data matrix to be as small as possible. The rationale is that we only expect a limited number of distinct cell types in the samples. The third term operates on the bulk RNA-Seq data. It ensures consistency between the average gene expression of the aggregated imputed data and the average gene expression of the bulk RNA-Seq data. We developed a convex optimization algorithm to minimize the objective function.

## R Version
### Install from Github 
```
library(devtools)
install_github("XXXXX")
```

### Install from source codes

* Download source codes [here](https://chopri.box.com/s/5vhh5ex450ml2pqf7ili5trzqnp5fumc) 
and In R type:
 
```
install.packages(path_to_file, type = 'source', rep = NULL)
```
Where path_to_file would represent the full path and file name:
- On Windows it will look something like this: "C:\\Download\SCRABBLE.tar.gz".
- On UNIX it will look like this: "~/Downloads/SCRABBLE.tar.gz".


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
## MATLAB Version

### Quick start
#### Load the data
There are three datasets in the .mat file. There are the true data set, Drop-out data set, and the imputed data set by SCRABBLE.
```
load('demo_data_HF.mat')
```
#### Prepare the data
We construct the data structure which is taken as one of the input of SCRABBLE.
```
data.data_sc = data_sc;
data.data_bulk = data_bulk;
```

#### Prepare the parameter for SCRABBLE
Set up the parameters used in example
```
parameter = [10,1e-5,1e-4];
nIter = 100;
```
#### Run SCRABBLE
```
dataRecovered = scrabble(data,parameter,nIter);
```
#### Plot the results
```
gcf = figure(1);
set(gcf, 'Position', [100, 500, 1200, 300])
subplot(1,3,1)
imagesc(log10(data_true+1))
title('True Data')
axis off
subplot(1,3,2)
imagesc(log10(data_sc+1))
title('Drop-out Data')
axis off
subplot(1,3,3)
imagesc(log10(dataRecovered+1))
title('Imputed Data by SCRABBLE')
axis off
```

## Help
Please feel free to contact Tao Peng (pengt@email.chop.edu) if you have any questions about the software.
## Reference
Tao Peng et al. Single-cell RNA-Seq imputation constrained by bulk RNA-Seq data. Submitted. 2018
