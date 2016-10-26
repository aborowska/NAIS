#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "blas.h"
#include "lapack.h"

#define  PI     3.14159265358979323846
#define Log2PI  1.83787706640934548356
#define LogPI   1.14472988584940017414

/* ********************************************************************** */

void my_gamma(double *gam, mwSignedIndex N)
{
    mxArray *in_array_ptr, *out_array_ptr; // mxArray * - a pointer to a struct (A POINTER TO A POINTER??)
   
    in_array_ptr = mxCreateDoubleMatrix(3*N, 1, mxREAL);  
    
    memcpy(mxGetPr(in_array_ptr), gam, 3*N*sizeof(double)); // start copying at the double returned by mxGetPr(array_ptr)
    mexCallMATLAB(1, &out_array_ptr, 1, &in_array_ptr, "gamma"); // & turns a value into a pointer --> call a MATLAB function with pointer to real matrices
    memcpy(gam, mxGetPr(out_array_ptr), 3*N*sizeof(double)); // start copying at the double returned by d=mxGetPr(plhs[0])
 
    mxDestroyArray(in_array_ptr);
    mxDestroyArray(out_array_ptr);
}

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

void EIS_reg_vec_copula(double *y, double *theta_smooth, double *V_smooth, double *z, double *w, double *tol_C,
                 double *nu, mwSignedIndex st, double *link,
                 mwSignedIndex T, mwSignedIndex M,
                 double *b_new, double *C_new)
{
    mwSignedIndex i, j;
    double *theta_GH, *Y, beta[2], gam[3], rho;
     	
    /* Variable size arrays */
    theta_GH = malloc((M)*sizeof(double));              
    Y = malloc((M)*sizeof(double));   
    
    if (st == 1)
    {        
        gam[0] = (nu[0]+2)/2;
        gam[1] = nu[0]/2;            
        gam[2] = (nu[0]+1)/2;    
        my_gamma(&gam,1);
//         mexPrintf("gam[0] = %6.4f\n",gam[0]);
//         mexPrintf("gam[1] = %6.4f\n",gam[1]);
//         mexPrintf("gam[2] = %6.4f\n",gam[2]);

    }
    
    for (i=0; i<T; i++)
    {
        for (j=0; j<M; j++)
        {
            theta_GH[j] = theta_smooth[i] + sqrt(V_smooth[i])*z[j];
            
            /* rho obtained from theta via link function */
            if (link[0]) /* KLS link*/
            {
                rho = (1 - exp(-theta_GH[j]))/(1 + exp(-theta_GH[j]));    
//                 mexPrintf("KLS link\n");
            }
            else    /* HM link*/
            {
                rho = (exp(2*theta_GH[j])-1)/(exp(2*theta_GH[j])+1);        
//                 mexPrintf("HM link\n");            
            }                  
                
            if (st == 1) /* Student's t distribution*/
            {
                Y[j] = y[i]*y[i] + y[i+T]*y[i+T] - 2*rho*y[i]*y[i+T];
                Y[j] = Y[j]/(nu[0]*(1-rho*rho));
                Y[j] = (nu[0] + 2)*log(1 + Y[j]);
                Y[j] = Y[j] + log(1 - rho*rho) - (nu[0] + 1)*(log(1 + y[i]*y[i]/nu[0]) + log(1 + y[i+T]*y[i+T]/nu[0]));  
                Y[j] = -0.5*Y[j];                   
                Y[j] = Y[j] + log(gam[0]);           
                Y[j] = Y[j] + log(gam[1]); 
                Y[j] = Y[j] - 2*log(gam[2]);                 
            }
            else
            {
                Y[j] = y[i]*y[i] + y[i+T]*y[i+T] - 2*rho*y[i]*y[i+T];
                Y[j] = Y[j]/(1-rho*rho);
                Y[j] = Y[j] + log(1-rho*rho) - y[i]*y[i] - y[i+T]*y[i+T]; 
                Y[j] = -0.5*Y[j]; 
            }      
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
    free(theta_GH); 
    free(Y);      
}

/* ********************************************************************** */


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex T, M, st;                                     /* size of matrix */
    double *y, *theta_smooth, *V_smooth, *z, *w, *tol_C;        /* input*/
    double *nu, *link;
    double *b_new, *C_new;                                      /* output */
    
    /* Getting the inputs */
    y = mxGetPr(prhs[0]);
    theta_smooth = mxGetPr(prhs[1]);
    V_smooth = mxGetPr(prhs[2]);
    z = mxGetPr(prhs[3]);
    w = mxGetPr(prhs[4]);
    tol_C = mxGetPr(prhs[5]);
	nu = mxGetPr(prhs[6]);
    link = mxGetPr(prhs[7]);
            
    T = mxGetM(prhs[0]); 
    M = mxGetM(prhs[3]); 
    st = mxGetM(prhs[6]);        

//     mexPrintf("st = %i\n", st);            
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(T,1,mxREAL);  
    plhs[1] = mxCreateDoubleMatrix(T,1,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    b_new = mxGetPr(plhs[0]);
    C_new = mxGetPr(plhs[1]);
   
    /* call the function */
    EIS_reg_vec_copula(y, theta_smooth, V_smooth, z, w, tol_C, nu, st, link,
            T, M, b_new, C_new);
  
}
