function [par_SV, hessian, hessian_tr, theta_smooth] = estimate_NAIS(par_SV_init, y, fn_jacobian, cont, options)
    n = length(y);
    cont.S = 200; % number of simulated trajectories in IS estimation
    RND = randn(n,S/2);   % normal random numbers; used in SimSmooth, for S/2 simulation paths; 
 
    %%  NAIS algorithm for selecting the importance parameter set [b,C]
    par_init.b = zeros(n,1);           
    par_init.C = ones(n,1);           % V = (Z*P1*Z')';
    par_NAIS = NAIS_param(par_init, y, par_SV_init, cont); % Algorithm 2: Efficient importance parameters via NAIS

    %% MAIN optimisation
    par_SV_init_trans = transform_param_SV(par_SV_init,  [cont.data_on,'_opt']);
    NAIS_max = @(par_SV_trans) NAIS_loglik(par_SV_trans, par_NAIS, y, S, cont, RND);
    
    [par_SV_trans,~,~,~,~, hessian]= fminunc(NAIS_max, par_SV_init_trans, options);
 
    jaco_inv = fn_jacobian(theta_trans);
    hessian_tr = jaco_inv*(n*hessian)*jaco_inv; 
    
    if (nargout == 4)
        [~, theta_smooth] = NAIS_max(par_SV_trans);
    end
    
    par_SV = transform_param_SV(par_SV_trans,  [cont.data_on,'_back']);

end