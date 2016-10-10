#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "blas.h"
#include "lapack.h"


void EIS_reg(double *y, double *theta, double *w, mwSignedIndex S,
           double *beta)
{
    mwSignedIndex i;
    char *chN = "N", *chT = "T";
    double one = 1.0, zero = 0.0;
    double A[9], B[3], *X, *Y;
    mwSignedIndex K, N;
    mwSignedIndex IPIV[3], INFO;

    /* Variable size arrays */
    X = malloc((3*S)*sizeof(double));              
    Y = malloc((S)*sizeof(double));              
    
    K = 3;
    N = 1;
    
    /* create X and Y */
    for (i=0; i<S; i++)
    {      
        X[i] = w[i];
        X[S+i] = w[i]*theta[i];
        X[2*S+i] = -0.5*theta[i]*X[S+i];
        Y[i] = w[i]*y[i];
    }    
    
    dgemm(chT, chN, &K, &K, &S, &one, X, &S, X, &S, &zero, &A, &K);      /* get X'*X*/
    dgemm(chT, chN, &K, &N, &S, &one, X, &S, Y, &S, &zero, &B, &K);      /* get X'*Y*/
    dgetrf(&K, &K, &A, &K, &IPIV, &INFO);                                /* get LU factorisation of A*/
    dgetrs(chN, &K, &N, &A, &K, &IPIV, &B, &K, &INFO);                   /* get solution B*/
    
    beta[0] = B[1];
    beta[1] = B[2];
    
    /* Free allocated memory */
    free(X); free(Y);    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex S;                                    /* size of matrix */
    double *y, *theta, *w;                              /* input*/
    double *beta;                                       /* output */
    
    /* Getting the inputs */
    y = mxGetPr(prhs[0]);
    theta = mxGetPr(prhs[1]);
    w = mxGetPr(prhs[2]);

    S = mxGetM(prhs[0]); 
             
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(2,1,mxREAL);  

    /* get a pointer to the real data in the output matrix */
    beta = mxGetPr(plhs[0]);
    
    /* call the function */
    EIS_reg(y, theta, w, S, beta);
  
}
