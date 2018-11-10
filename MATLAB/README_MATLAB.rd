SCRABBLE is a scRNA-Seq imputation method by combining scRNA-Seq and 
bulk RNA-Seq data. 

***************************************************************************
***************************** INTRODUCTION ********************************
***************************************************************************

Single cell RNA-Seq (scRNA-Seq) data suffers from a large proportion of 
zeros or low read counts for expressed genes. Such drop-out events present 
a fundamental challenge for downstream analysis. Here we describe the 
SCRABBLE (Single Cell RNA-Seq imputAtion constrained By BuLk RNAsEq data) 
algorithm to address this issue. We demonstrate that SCRABBLE outperforms 
existing methods in recovering drop-out events and preserving information 
about gene-gene relationship and cell-cell relationship in the data.

***************************************************************************
*************************** MATHEMATICAL MODEL ****************************
***************************************************************************

X^* = argmin_{X >=0}(1/2||P_\Omega(X) - X_0||_F^2 + \alphaRank(X) + \beta||AX - D||)

where X_0 is the drop-out scRNAseq data set. D is the bulk RNAseq data. 
\Omega and A are known. \alpha and \beta are two parameters in the model. 

(The details of the mathematics deduction are shown in the manuscript)

***************************************************************************
********************************* INPUT ***********************************
***************************************************************************

Prepare data: 
The purpose of the algorithm is used to impute the scRNAseq data. The common 
formats of scRNAseq are calculated as TPM or counts of UMI. Both of them 
could be used as the input of SCRABBLE.  

data: data is a structure data type. It contains two types of data,
scRNA-Seq and bulk RNA-Seq data, which are consistently normalized data.
The data matrix of scRNAseq data has the gene as the row and the cell as 
the column.
The bulk RNA-seq data is the column vector of the average expression of 
genes across all the bulk samples.
The number of rows of the scRNAseq and bulk RNAseq is the same.
To reduce the computation cost, we recommend to reduce the number of genes 
in the data using the experesion levels of the genes as the standard. 
In the real data, we could get around 10,000 genes after setting the 
threshold as 10 UMIs. 

parameter: parameter is a vector containing two elements. The first element
is the alpha parameter, which is the weight for matrix rank. 
The second element is the beta parameter, which is the weight for 
the consistency between the aggregated scRNA-Seq data and 
the bulk RNA-Seq data. The recommendation range of the \alpha is from 1 
to 500. The recommendation range of \beta is proportional to \alpha and 
the size of the input data matrix. 

nIter: nIter is the maximal number of iterations. The recommmdation times 
of iteration are 100.

***************************************************************************
******************************** OUTPUT ***********************************
***************************************************************************

The output is the imputed data matrix. Each row stands one gene. 
Each column is one cell. The orders of genes and cells are the same with 
original data matrix. It could be ready for the downstream analysis as new 
TPM and UMIs directly. 

***************************************************************************
******************************* EXAMPLE ***********************************
***************************************************************************

"demo_simulation.m" is an example to show how to use the SCRABBLE software. 
You could input "demo_simulation" at MATLAB command window to run the example.
"html" folder contains the html file which shows the demo result.  

***************************************************************************
****************************** CONTACT ************************************
***************************************************************************

Please feel free to contact Tao Peng (pengt@email.chop.edu) if you have any 
questions about the software.
