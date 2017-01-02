/*
 * function [x1,x2]=prox_lse_2D(a,y1,y2)
 *
 * solves
 * argmin_(x1,x2) a*log(exp(x1)+exp(x2)) + 1/2(x1-y1)^2 + 1/2(x2-y2)^2
 *
 * INPUT :
 * - a : real >0
 * - y1,y2 : arrays of real numbers.
 *
 * OUTPUT :
 * - x1,x2 : arrays of the same size as y1,y2.
 *
 * Developer : P. Weiss December 2016.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "mex.h"
#include <omp.h>

void prox_lse_2D(double *a,int n, double *y1, double *y2, double* x1,double* x2){
    double tol1=1e-16;
    double tol2=1e-15;
    int nitmax=100;
    
    int unsolved;
    double z1,z2,lambda,aa;
    double f,fp,d;
    int i,k;
    
    #pragma omp parallel for private(z1,z2,f,fp,d,i,k,aa,unsolved,lambda)
    for (i=0;i<n;++i){
        
        z1=y1[i];z2=y2[i];
        aa=a[i];
        unsolved=1;
        
        if (z1>=z2){
            
            if (z2-z1+aa<=log(tol1)){
                lambda=1;
                unsolved=0;
            } else{
                lambda=1.0/(1+exp(z2-z1))-tol1;
            }
            
            k=0;
            while (unsolved && k<nitmax){
                f=z2-z1-aa+2*aa*lambda+log(lambda/(1.0-lambda));
                fp=2*aa+1.0/lambda+1.0/(1.0-lambda);
                d=f/fp;
                lambda=lambda-d;
                ++k;
               
                if (fabs(d)<tol2){unsolved=0;}
            }
            
            x1[i]=z1-aa*lambda;x2[i]=z2-aa*(1-lambda);
        } else{
            
            if (z1-z2+aa<=log(tol1)){
                lambda=1;
                unsolved=0;
            } else{
                lambda=1.0/(1.0+exp(z1-z2))-tol1;
            }
            
            k=0;
            while (unsolved && k<nitmax){
                f=z1-z2-aa+2*aa*lambda+log(lambda/(1.0-lambda));
                fp=2*aa+1.0/lambda+1.0/(1.0-lambda);
                d=f/fp;
                lambda=lambda-d;               
                ++k;
                
                if (fabs(d)<tol2){unsolved=0;}
            }
            x1[i]=z1-aa*(1-lambda);x2[i]=z2-aa*lambda;
        }
    }
}

/***************************************************************************/
/*En-tete Matlab du fichier                                                 */
/* Arguments a donner dans l'ordre a,z1,z2 */
/* Output : x1,x2 */
/* compile with  mex prox_lse_2D.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" */
/* Launch Matlab with OMP_NUM_THREADS=8 matlab& */
/* Or equivalently type maxNumCompThreads(N) before calling the function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    /*Declarations*/
    int nelements;
    double *y1,*y2,*a;
    double *x1,*x2;
    
    // Set up openmp
    int nProcessors=omp_get_max_threads();
    //mexPrintf("Number of processors: %d\n",nProcessors);
    omp_set_num_threads(nProcessors);
    
    
    /*Check for proper number of input arguments*/
    if(nrhs!=3) {
        mexErrMsgTxt("3  input required.");
    } else if(nlhs>3) {
        mexErrMsgTxt("Too many  output arguments");
    }
    
    /*Get the number of elements in the input argument*/
    nelements=mxGetNumberOfElements(prhs[1]);
    
    /*Get the input argument*/
    a =mxGetPr(prhs[0]);
    y1=mxGetPr(prhs[1]);
    y2=mxGetPr(prhs[2]);
    
    /*Allocates output arguments*/
    plhs[0]=mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]),
            mxGetDimensions(prhs[1]),
            mxGetClassID(prhs[1]),
            mxREAL);
    
    plhs[1]=mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]),
            mxGetDimensions(prhs[1]),
            mxGetClassID(prhs[1]),
            mxREAL);
    
    x1 = mxGetPr(plhs[0]);
    x2 = mxGetPr(plhs[1]);
    
    /* Main function */
    prox_lse_2D(a,nelements,y1,y2,x1,x2);
}
