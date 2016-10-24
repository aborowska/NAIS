% This is a very early ilustration of the NAIS algorithm (a really explicit one). 
% A more up-to-date version is in estimate_SV.m file (based on MEXed functions).
% No parameter estimation, just construction of the optimal IS parameters
% given the chosen SV model parameters.
% There is also some plots generation.

clear all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

par_SV.c = 1;
par_SV.phi = 0.98;
par_SV.sigma2 = 0.0225;
par_SV.n = 2000;

y = sim_SV(par_SV); % simulated daily log-returns
figure(1)
set(gcf,'defaulttextinterpreter','latex');
plot(y)
plotTickLatex2D

n = length(y);
Un = ones(n,1); 

%%
% NAIS algorithm for selecting the importance parameter set [b,C]
par_init.b = zeros(n,1);
par_init.C = ones(n,1);
cont.M = 30; % number of the Gauss-Hermite nodes
cont.tol = 0.0001; % convergence tolerance
cont.fast = true;

% Algorithm 2: Efficient importance parameters via NAIS
par_NAIS = NAIS_param(cont,par_init,y,par_SV);
b = par_NAIS.b;
C = par_NAIS.C;
a = 0.5*(log(abs(C)) - 2*log(pi) - (b.^2)./C);
y_star = b./C;

figure(2)
set(gcf,'defaulttextinterpreter','latex');
subplot(2,1,1)
plot(b)
title('Importance parameter b')
plotTickLatex2D
subplot(2,1,2)
plot(C)
title('Importance parameter C')
plotTickLatex2D

%%
% Algorithm 1: IS using an approximation linear state model
% (with antithetic variables for variance reduction)

% Given the optimal IS parameters obtain the smoothed mean of signal for
% the importance model 
par_KFS = IS_Model(par_NAIS, par_SV);
[theta_smooth, V, nu, F] = KFS(y_star,par_KFS); 

figure(3)
% set(gcf,'defaulttextinterpreter','latex');
plot(y)
hold on
plot(theta_smooth,'r')
hold off
% title('Observation $$y_{t}$$ and smoothed signal $$\hat{\theta}_t$$');
legend({'$$y_{t}$$','$$\hat{\theta}_t$$'},'Location','SouthEast','interpreter', 'latex');
plotTickLatex2D



%% Importance sampling:
% simulation smoothing for linear state space models to sample S/2
% trajectories for the signal e.g. via JSDK
S = 500;
theta_sim = zeros(n,S);
theta_sim(:,1:S/2) = IS_sim(S,n,theta_smooth,par_KFS);
% antithetic draws
theta_sim(:,(S/2+1):S) = 2*repmat(theta_smooth,1,S/2) - theta_sim(:,1:S/2);

%% LogLikelihood evaluation
% loglikelihood of the approximation linear state space model via e.g. KF
lng_y = -0.5*(n*log(2*pi) + sum(log(F)) + sum((nu.^2)./F));

% compute the logweights
lnP =  -0.5*(log(2*pi) + log(exp(par_SV.c + theta_sim)) - repmat(y.^2,1,S).*exp(-(par_SV.c + theta_sim)));
lnG = repmat(a,1,S) + repmat(b,1,S).*theta_sim -0.5*repmat(C,1,S).*theta_sim.^2;
logp = sum(lnP, 1);
logg = sum(lnG, 1);
lnw_s = logp - logg; 
lnw_s = lnw_s - max(lnw_s); % robustification
w_s = exp(lnw_s);
% compute the loglikelihood estimate
lnL_hat = lng_y + log(sum(w_s)) - log(S); 
fprintf('The estimated loglikelihood is: %6.4f.',lnL_hat)
