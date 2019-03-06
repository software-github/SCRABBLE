
[![DOI](https://zenodo.org/badge/166503417.svg)](https://zenodo.org/badge/latestdoi/166503417)

# SCRABBLE
<b>S</b>ingle <b>C</b>ell <b>R</b>NA-Seq imput<b>A</b>tion constrained <b>B</b>y <b>B</b>u<b>L</b>k RNAs<b>E</b>q data (SCRABBLE)

## SCRABBLE has been implemented in R and MATLAB.

SCRABBLE imputes drop-out data by optimizing an objective function that consists of three terms. The first term ensures that imputed values for genes with nonzero expression remain as close to their original values as possible, thus minimizing unwanted bias towards expressed genes. The second term ensures the rank of the imputed data matrix to be as small as possible. The rationale is that we only expect a limited number of distinct cell types in the samples. The third term operates on the bulk RNA-Seq data. It ensures consistency between the average gene expression of the aggregated imputed data and the average gene expression of the bulk RNA-Seq data. We developed a convex optimization algorithm to minimize the objective function.

## R Version
### Install from Github 
```
library(devtools)
install_github("software-github/SCRABBLE/R")
```

### Install from source codes

Download source codes [here](https://github.com/software-github/SCRABBLE/blob/master/SCRABBLE_0.0.1.tar.gz?raw=true) 


and In R type:
 
```
install.packages(path_to_file, type = 'source', rep = NULL)
```
Where `path_to_file` would represent the full path and file name:
- On Windows it will look something like this: "C:\\Downloads\SCRABBLE.tar.gz".
- On UNIX it will look like this: "~/Downloads/SCRABBLE.tar.gz".


### Quick start
```
data_sc <- demo_data[[1]]
data_bulk <- demo_data[[2]]
data_true <- demo_data[[3]]

parameter <- c(1,1e-6,1e-4)

result <- scrabble(demo_data, parameter = parameter)
```                   
## MATLAB Version

### Quick start
#### Load the data
There are three datasets in the .mat file. There are the true data set, Drop-out data set, and the imputed data set by SCRABBLE.
```
load('demo_data.mat')
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
parameter = [1,1e-6,1e-4];
```
#### Run SCRABBLE
```
dataRecovered = scrabble(data,parameter);
```
#### Visualize the results
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
Please feel free to contact Tao Peng (software.github@gmail.com) if you have any questions about the software.
## Reference
Tao Peng et al. Single-cell RNA-Seq imputation constrained by bulk RNA-Seq data. Submitted. 2018

