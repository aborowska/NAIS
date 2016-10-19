#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "blas.h"
#include "lapack.h"


/* ********************************************************************** */

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

/* ********************************************************************** */

void EIS_reg_vec_t(double *y, double *theta_smooth, double *V_smooth, double *z, double *w, double *tol_C,
                 double *nu, double *pdf_const,
				 mwSignedIndex T, mwSignedIndex M,
                 double *b_new, double *C_new)
{
    mwSignedIndex i, j;
    double *theta_GH, *Y, beta[2];
     	
    /* Variable size arrays */
    theta_GH = malloc((M)*sizeof(double));              
    Y = malloc((M)*sizeof(double));              

    for (i=0; i<T; i++)
    {
        for (j=0; j<M; j++)
        {
            theta_GH[j] = theta_smooth[i] + sqrt(V_smooth[i])*z[j];
            Y[j] = pdf_const[0] - 0.5*(theta_GH[j] + (nu[0]+1)*log(1 + y[i]*y[i]/((nu[0]-2)*exp(theta_GH[j]))));
        }
        
        EIS_reg(Y, theta_GH, w, M, beta);   
        
        b_new[i] = beta[0];
        if (beta[1] < *tol_C)
        {
            C_new[i] = *tol_C;
        }
        else
        {
            C_new[i] = beta[1];
        }
    
    }

    /* Free allocated memory */
    free(theta_GH); free(Y);  

}

/* ********************************************************************** */


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex T, M;                                            /* size of matrix */
    double *y, *theta_smooth, *V_smooth, *z, *w, *tol_C;        /* input*/
    double *nu, *pdf_const;
    double *b_new, *C_new;                                       /* output */
    
    /* Getting the inputs */
    y = mxGetPr(prhs[0]);
    theta_smooth = mxGetPr(prhs[1]);
    V_smooth = mxGetPr(prhs[2]);
    z = mxGetPr(prhs[3]);
    w = mxGetPr(prhs[4]);
    tol_C = mxGetPr(prhs[5]);
	nu = mxGetPr(prhs[6]);
 	pdf_const = mxGetPr(prhs[7]);
    
    T = mxGetM(prhs[0]); 
    M = mxGetM(prhs[3]); 
             
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(T,1,mxREAL);  
    plhs[1] = mxCreateDoubleMatrix(T,1,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    b_new = mxGetPr(plhs[0]);
    C_new = mxGetPr(plhs[1]);
    
    /* call the function */
    EIS_reg_vec_t(y, theta_smooth, V_smooth, z, w, tol_C, nu, pdf_const, T, M, b_new, C_new);
  
}
