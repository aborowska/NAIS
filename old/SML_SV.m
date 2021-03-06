clear all
addpath(genpath('../include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

data_on = true;
% data_on = false;


%% NAIS control parameters
cont.M = 20; % number of the Gauss-Hermite nodes
cont.tol = 0.0001; % convergence tolerance
cont.tol_C = 1e-5; % tolerance for the variance 
[cont.GH.z, cont.GH.h]  = hernodes(cont.M);
cont.GH.z = -cont.GH.z;
cont.GH.h = cont.GH.h.*exp(0.5*(cont.GH.z).^2);
cont.iter_max = 20;
% cont.err = 'n';
cont.S = 200;  % number of simulated trajectories in IS estimation
% then:  RND are normal random numbers, used in SimSmooth, for S/2 simulation paths; 
% 2 columns for each simulation; S/2 simulations are run one column for eta and one for epsilon; 
if data_on
    cont.data_on = 'est'; 
else
    cont.data_on = 'sim';
end 
cont.print = false; % true for printing iteration info in estimation  


%% Observations
if data_on
    % Use data
%     y = csvread('GSPC_ret.csv');
%     y = 100*y;
% ind_arch = find(y<=-5.5, 1, 'last' );
% y = y(1:ind_arch,1);
% y = y - mean(y);

% %     y = csvread('IBM_ret.csv');
    y = csvread('GSPC_ret_updated_short_sv.csv'); % the most recent 5 years
%     y = csvread('GSPC_ret_updated.csv'); % crisis data    
    
    y = 100*y;
    n = length(y);
else    
    % Simulation
%     par_SV_sim = [0.5, 0.98, 0.15^2];
    par_SV_sim = [-10, 0.95, 0.25];
%     n = 1000;
%     n = 2500;
    n = 5000;
    y = sim_SV(par_SV_sim,n); % simulated daily log-returns
    
cont.data_on = 'est'; 
    
end

% MC = 20;
% par_SV_MC = zeros(MC,3);
% 
% for ii = 1:MC
    %% Optimisation
%     S = 200; % number of simulated trajectories in IS estimation
%     RND = randn(n,S/2);   % normal random numbers; used in SimSmooth, for S/2 simulation paths; 
                        % 2 columns for each simulation; S/2 simulations are run
                        % one column for eta and one for epsilon; 
    options=optimset('display','iter','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'maxiter',500,'HessUpdate','bfgs','FinDiffType','central');

    % Initial optimisation
    par_SV_init = [  0.9080    0.8910    0.0044];
%     par_SV_init = [0.5, 0.98, 0.15^2];
% par_SV_init = par_SV_sim;
%     par_SV_init = [0.1, 0.97, 0.03];
  
    try 
        [par_SV_adapt, hess_SV_adapt, hess_SV_corr_adapt] = NAIS(par_SV_init, y, S, cont, RND, options);
    catch 
        par_SV_adapt = par_SV_init;
    end
    % Final optimisation
    try
        [par_SV_opt, hess_SV_opt, hess_SV_corr_opt, theta_smooth] = NAIS(par_SV_adapt, y, S, cont, RND, options);
    catch
        par_SV_opt = par_SV_init;
    end

%     par_SV_MC(ii,:) = par_SV_opt;
% end

V_SV_opt = inv(hess_SV_opt);
std = sqrt(diag(V_SV_opt));

V_SV_corr_opt = inv(hess_SV_corr_opt);
std_corr = sqrt(diag(V_SV_corr_opt));


% SML_barra.par_SV_opt = par_SV_opt;
% SML_barra.V_SV_corr_opt = V_SV_corr_opt;
% 
% SML_tutorial.par_SV_opt = par_SV_opt;
% SML_tutorial.V_SV_corr_opt = V_SV_corr_opt;
% 
% SML_nais.par_SV_opt = par_SV_opt;
% SML_nais.V_SV_corr_opt = V_SV_corr_opt;

if data_on
% %     save 'results/SML_arch.mat' 'par_SV_opt' 'V_SV_corr_opt';    
% %     save 'results/SML_ibm.mat' 'par_SV_opt' 'V_SV_corr_opt';  
   save 'results/SML_gspc_updated_short_sv.mat' 'par_SV_opt' 'V_SV_corr_opt' 'theta_smooth';  
%    save 'results/SML_gspc_updated.mat' 'par_SV_opt' 'V_SV_corr_opt' 'theta_smooth';  
else
    save 'results/SML_sim.mat' 'par_SV_opt' 'V_SV_corr_opt';
end