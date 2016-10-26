clear all
addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

data_on = true; % data_on = false;
model = 'sv'; % 'svt'

%% NAIS control parameters
cont = NAIS_control(data_on);

%% Observations
if data_on % Use data
%   y = csvread('IBM_ret.csv');
    y = csvread('GSPC_ret_updated_short_sv.csv'); % the most recent 5 years
%   y = csvread('GSPC_ret_updated.csv'); % crisis data    
    y = 100*y;
    n = length(y);
else % Simulation
    if strcmp(model, 'sv')
%         par_SV_sim = [0.5, 0.98, 0.15^2];
        par_SV_sim = [-10, 0.95, 0.25]; 
    else
        par_SV_sim = [0.5, 0.98, 0.15^2, 10];
    end
    n = 5000;
    y = sim_volatility(par_SV_sim,n); % simulated daily log-returns
end

%% Initialisation 
options=optimset('display','iter','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'maxiter',500,'HessUpdate','bfgs','FinDiffType','central');

%     par_SV = [c, phi, sigma2_eta, nu]
if strcmp(model,'sv')
    par_SV_init = [  0.9080    0.8910    0.0044];
%         par_SV_init = [0.5, 0.98, 0.15^2];
%         par_SV_init = par_SV_sim;
%         par_SV_init = [0.1, 0.97, 0.03];
else
    par_SV_init = [0.5, 0.98, 0.15^2, 10];
end

fn_jacobian = @(xx) jacobian_ss(xx); % Jacobian of the parameter tranformation to get standard errors of the orignal parameters
    
%% Optimisation    
% Uncommnet the loop if MC replications of SML required
% MC = 20;
% par_SV_MC = zeros(MC,3);
% 
% for ii = 1:MC
    % Initial optimisation
    try 
        [par_SV_adapt, hess_SV_adapt, hess_SV_corr_adapt] = estimate_NAIS(par_SV_init, y, fn_jacobian, cont, options);
    catch 
        par_SV_adapt = par_SV_init;
    end
    % Final optimisation
    try
        [par_SV_opt, hess_SV_opt, hess_SV_corr_opt, theta_smooth] = estimate_NAIS(par_SV_adapt, y, fn_jacobian, cont, options);
    catch
        par_SV_opt = par_SV_init;
    end

%     par_SV_MC(ii,:) = par_SV_opt;
% end

V_SV_opt = inv(hess_SV_opt);
std = sqrt(diag(V_SV_opt));

V_SV_corr_opt = inv(hess_SV_corr_opt);
std_corr = sqrt(diag(V_SV_corr_opt));

% Choose the correct file name! (corresponding to the dataset)
if data_on
%     save 'results/SML_ibm.mat' 'par_SV_opt' 'V_SV_corr_opt';  
   save 'results/SML_gspc_updated_short_sv.mat' 'par_SV_opt' 'V_SV_corr_opt' 'theta_smooth';  
%    save 'results/SML_gspc_updated.mat' 'par_SV_opt' 'V_SV_corr_opt' 'theta_smooth';  
else
    save 'results/SML_sim.mat' 'par_SV_opt' 'V_SV_corr_opt';
end