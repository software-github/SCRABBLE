/*==========================================================
 * cDescent.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier)
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		newX1 = cDescent(gamma,Y,B,Lambda,A,zones,newX);
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 * /*==========================================================

 *
 *
 * /* The computational routine */
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#define gamma_IN prhs[0]
#define y_IN prhs[1]
#define b_IN prhs[2]
#define lambda_IN prhs[3]
#define a_IN prhs[4]
#define zones_IN prhs[5]
#define newx_IN prhs[6]
#define x_OUT plhs[0]

void cDescent(double gamma, double *y, double *b, double *lambda, double *a, double *zones, double *newx, double *x, int m1, int n1)
{   /* allocate the variables */
    int i,j,k,index;
    double tmp;
    /* initalize the outpouts */
    for (i = 0; i < m1; i++){
        for(j = 0; j < n1; j++) {
            index = i + m1*j;
            x[index] = newx[index];
        }
    }    
    /* Outputs updated */
    for (i = 0; i < m1; i++){
        for(j = 0; j < n1; j++) {
            tmp = 0;
            for (k = 0; k < m1; k++) {
                tmp += a[i + m1*k]*x[k + m1*j];
            }
            index = i + m1*j;
            tmp = (gamma*y[index] + b[index] - lambda[index] - tmp + a[i + m1*i]*x[index])/(a[i + m1*i] + zones[index]);
            x[index] = max(tmp,0);
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double gamma;              /* input scalar */
    int m1,n1;                 /* intermeidate variable */
    double *y, *b, *lambda, *a, *zones, *newx;     /* input matrix */     
    double *x;              /* output matrix */
    
    /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) ||
            mxIsComplex(prhs[0]) ||
            mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) ||
            mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[2]) ||
            mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[3]) ||
            mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[4]) ||
            mxIsComplex(prhs[4])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[5]) ||
            mxIsComplex(prhs[5])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[6]) ||
            mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }    
    
    /* get the value of the scalar input  */
    gamma = mxGetScalar(gamma_IN);
    /* get the value of the double input matrix*/
    m1 = mxGetM(y_IN);  /* Get pointer to Y's row number */
    n1 = mxGetN(y_IN);  /* Get pointer to Y's column number */
    /* Assign pointers to the intput matrix */
    y = mxGetPr(y_IN);  /* Get pointer to Y's data */
    a = mxGetPr(a_IN);  /* Get pointer to A's data */
    b = mxGetPr(b_IN);  /* Get pointer to B's data */
    lambda = mxGetPr(lambda_IN);  /* Get pointer to Lambda's data */
    zones = mxGetPr(zones_IN);    /* Get pointer to zones's data */
    newx = mxGetPr(newx_IN);      /* Get pointer to newX's data */
    
    /* Create matrix for the return argument. */
    x_OUT = mxCreateDoubleMatrix(m1, n1, mxREAL);
    /* Assign pointers to the output matrix */ 
    x = mxGetPr(x_OUT);
    
    /* Do the actual computations in a subroutine */ 
    cDescent(gamma,y,b,lambda,a,zones,newx,x,m1,n1);
}
