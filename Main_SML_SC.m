clear all
addpath(genpath('include/'));

%     location = 'C:\Program Files\MATLAB\R2014a\extern\lib\win32\lcc\';
% %     location = 'C:\Program Files\MATLAB\R2014a\extern\lib\win64\microsoft\';
%     lapacklib = [location 'libmwlapack.lib'];
%     blaslib = [location 'libmwblas.lib'];
%     mex('-v', '-largeArrayDims', 'EIS_reg_vec_copula.c', blaslib, lapacklib) 
      
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

data_on = false; % data_on = true;
model = 'sct'; % 'sct'

%% NAIS control parameters
cont = NAIS_control(data_on);
cont.link = 1;  % 1: 2*(logsig-0.5) (KLS); 0: inv Fisher transform (Hafner & Manner 2012);
cont.print = false;

%% Observations
if data_on % Use data
    y = load('ibm_ccola_rets.txt'); % for Patton's copula toolbox
    T = size(y,2);
else % Simulation
    if strcmp(model, 'sc')
        par_sim = [1, 0.98, 0.01];       
    else
        par_sim = [1, 0.98, 0.01, 10];       
    end
    T = 2500;
    N = 1;
    link = 1;
    [y, theta, rho] = sim_copula(par_sim, T, N, link);
end

%% Initialisation 
%     par_SV = [c, phi, sigma2_eta, nu]
%         par_SC_init = par_sim;
if strcmp(model,'sc')
    par_SC_init = [0.9, 0.97, 0.15^2]; 
    par_SC_init = [1, 0.98, 0.01];
else
    par_SC_init = [1, 0.97, 0.15^2, 10];

    par_SC_init = [1, 0.98, 0.01, 10];
end

options = optimset('display','iter','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'maxiter',500,'HessUpdate','bfgs','FinDiffType','central');
fn_jacobian = @(xx) jacobian_ss(xx); % Jacobian of the parameter tranformation to get standard errors of the orignal parameters
RND = randn(T,cont.S/2);   % normal random numbers; used in SimSmooth, for S/2 simulation paths; 

%% Optimisation    
% Uncommnet the loop if MC replications of SML required
% MC = 20;
% par_SC_MC = zeros(MC,3);
% 
% for ii = 1:MC

% KLS = load('kls_sim_copula.mat');       
%     y=y_kls(:,1:2500)';
GAS = load('Results/ibm_ccola_copula_gas_results.mat', 'u', ...
    'mu_copula', 'Sigma_copula', 'f_copula', 'rho_copula',...
    'mu_copula_t', 'Sigma_copula_t', 'f_copula_t', 'rho_copula_t');
u = GAS.u;

    % Initial optimisation
    try 
        [par_SC_adapt, hess_SC_adapt, hess_SC_corr_adapt] = estimate_NAIS_copula(par_SC_init, y, fn_jacobian, cont, RND, options);
%           [par_SC_adapt, hess_SC_adapt, hess_SC_corr_adapt] = estimate_NAIS_copula(par_SC_init, u, fn_jacobian, cont, options);
    catch 
        par_SC_adapt = par_SC_init;
    end
    % Final optimisation
    try         
        [par_SC_opt, hess_SC_opt, hess_SC_corr_opt, theta_smooth] = estimate_NAIS_copula(par_SC_adapt, y, fn_jacobian, cont, RND, options);
%         [par_SC_opt, hess_SC_opt, hess_SC_corr_opt, theta_smooth] = estimate_NAIS_copula(par_SC_adapt, u, fn_jacobian, cont, options);
%par_SC_opt =   0.3934    0.9966    0.0012
%  KLS:  [0.3922, 0.9962, 0.0013]       
    catch
        par_SC_opt = par_SC_init;
    end

%     par_SC_MC(ii,:) = par_SC_opt;
% end

V_SC_opt = inv(hess_SC_opt);
std = sqrt(diag(V_SC_opt));

V_SC_corr_opt = inv(hess_SC_corr_opt);
std_corr = sqrt(diag(V_SC_corr_opt));

% Choose the correct file name! (corresponding to the dataset)
if data_on
    save 'results/SML_ibm_ccola_copula.mat' 'par_SC_opt' 'V_SC_corr_opt' 'theta_smooth';  
 else
    save 'results/SML_copula_sim.mat' 'par_SC_opt' 'V_SC_corr_opt';
end