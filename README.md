# NAIS
My implementation of Numerically Accelerated Importance Sampling (NAIS) for nonlinear non-Gaussian state space models by Koopman et al. (2015). NAIS is an efficient modification of Efficient Importance Sampling by Richard and Zhang (2007) based on numerical integration using the Gauss-Hermite quadrature instead of sampling based MC integration. The key assumption in NAIS is linear Gaussian state equation. 

I use it mainly for two purposes:
* SML (Simulated Maximum Likelihood) estimation for Stoachstic Volatility (SV) and SVt models;
* states sampling (given the model parameters) in Bayesian analysis.

To speed up computations, I "MEXed" some of the MATLAB codes (i.e. they are written in C).

### References
Koopman, S. J., A. Lucas and M. Scharth (2015), "Numerically Accelerated Importance Sampling for Nonlinear Non-Gaussian State Space Models", _Journal of Business and Economic Statistics_, 33, 114-127.

Richard, J. and W. Zhang (2007), "Efficient High-Dimensional Importance Sampling", __Journal of Econometrics__, 141, 1385-1411.
