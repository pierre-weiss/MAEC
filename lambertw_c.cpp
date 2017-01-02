/*
 * function w=lambertw(z)
 *
 * solves
 * w exp(w) = z
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

void lambertw(double *z,int nelements,double *w){

    int nitmax=36;    
    double c1,c2,w1,dw,tmp,wt;
    int i,k,solved;
    
    #pragma omp parallel for private(c1,c2,w1,dw,tmp,wt,k,solved,nitmax)
    for (i=0;i<nelements;++i){

        if (fabs(z[i] + 0.3678794411714423216) <= 1.5){
            wt = sqrt(5.43656365691809047*z[i] + 2) - 1;
        } else {
            tmp=log(z[i]+(z[i]==0));
            wt=tmp - log(tmp+(tmp==0));
        }
        
        
        k=1;solved=0;
        while(k<=nitmax && (solved==0)){
            c1=exp(wt);
            c2=wt*c1-z[i];
            w1=wt+(wt!=-1);            
            dw=c2/(c1*w1 - ((wt + 2)*c2/(2*w1)));
            wt=wt-dw;
            if (fabs(dw)<0.7e-16*(2+abs(wt))){
                solved=1;
            }
            ++k;
        }
        w[i]=wt;
    }
}

/***************************************************************************/
/*En-tete Matlab du fichier                                                 */
/* Input : z */
/* Output : w */
/* compile with  mex lambertw_c.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" */
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
    lambertw(z,nelements,w);
}
