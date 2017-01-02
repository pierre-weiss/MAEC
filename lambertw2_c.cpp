/*
 * function w=lambertw2(z)
 *
 * solves
 * w exp(w) = exp(z)
 *
 * INPUT :
 * - z: array of reals.
 *
 * OUTPUT :
 * - w: lambert of z.
 *
 * Developer : P. Weiss December 2016.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "mex.h"
#include <omp.h>

void lambertw2(double *z,int nelements,double *w){
    
    int nitmax=10;
    double tol=1e-15;
    double c1,c2,w1,dw,tmp,wt,t;
    int i,k,solved;
    
    #pragma omp parallel for private(c1,c2,w1,dw,tmp,wt,t,k,solved)
    for (i=0;i<nelements;++i){
        if (z[i]==-INFINITY){
            w[i]=0;
        } else if (z[i]<=0.12){
            t=exp(z[i]);
            wt = sqrt(5.43*t + 2.0) - 1.0;
            k=1;solved=0;
            while(k<=nitmax && (solved==0)){
                c1=exp(wt);
                c2=wt*c1-t;
                w1=wt+ (wt!=-1);
                dw=c2/(c1*w1 - ((wt + 2.0)*c2/(2*w1)));
                wt=wt-dw;

                if (fabs(dw)<tol*(2+fabs(wt))){
                    solved=1;
                }
                ++k;
            }
            w[i]=wt;
        } else {
            wt=z[i]-log(z[i]);
            k=1;solved=0;
            while(k<=nitmax && (solved==0)){
                c1=log(wt)+wt-z[i];
                c2=1.0/wt+1;
                dw=c1/c2;
                wt=wt-dw;
                if (fabs(dw)<tol*(2+fabs(wt))){
                    solved=1;
                }
                ++k;
            }
            w[i]=wt;
        }
    }
}

/***************************************************************************/
/*En-tete Matlab du fichier                                                 */
/* Input : z */
/* Output : w */
/* compile with  mex lambertw2_c.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" */
/* Launch Matlab with OMP_NUM_THREADS=8 matlab& */
/* Or equivalently type maxNumCompThreads(N) before calling the function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    /*Declarations*/
    int nelements;
    double *z;
    double *w;
    
    // Set up openmp
    int nProcessors=omp_get_max_threads();
    //mexPrintf("Number of processors: %d\n",nProcessors);
    omp_set_num_threads(nProcessors);
    
    
    /*Check for proper number of input arguments*/
    if(nrhs!=1) {
        mexErrMsgTxt("1  input required.");
    } else if(nlhs>1) {
        mexErrMsgTxt("Too many  output arguments");
    }
    
    /*Get the number of elements in the input argument*/
    nelements=mxGetNumberOfElements(prhs[0]);
    
    /*Get the input argument*/
    z=mxGetPr(prhs[0]);
    
    /*Allocates output arguments*/
    plhs[0]=mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
            mxGetDimensions(prhs[0]),
            mxGetClassID(prhs[0]),
            mxREAL);
    
    w = mxGetPr(plhs[0]);
    
    /* Main function */
    lambertw2(z,nelements,w);
}
