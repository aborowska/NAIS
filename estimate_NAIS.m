function [par_SV, hessian, hessian_tr, theta_smooth] = estimate_NAIS(par_SV_init, y, fn_jacobian, cont, options)
    n = length(y);
    RND = randn(n,cont.S/2);   % normal random numbers; used in SimSmooth, for S/2 simulation paths; 
 
    %%  NAIS algorithm for selecting the importance parameter set [b,C]
    par_NAIS_init.b = zeros(n,1);           
    par_NAIS_init.C = ones(n,1);           % V = (Z*P1*Z')';
    par_NAIS = NAIS_param(par_NAIS_init, y, par_SV_init, cont); % Algorithm 2: Efficient importance parameters via NAIS

    %% MAIN optimisation
    par_SV_init_trans = transform_param_ss(par_SV_init,  [cont.data_on,'_opt'], cont);
    NAIS_max = @(xx) NAIS_loglik(xx, par_NAIS, y, cont.S, cont, RND);
    
    [par_SV_trans,~,~,~,~, hessian]= fminunc(NAIS_max, par_SV_init_trans, options);
 
    jaco_inv = fn_jacobian(par_SV_trans);
    hessian_tr = jaco_inv*(n*hessian)*jaco_inv; 
    
    if (nargout == 4)
        [~, theta_smooth] = NAIS_max(par_SV_trans);
    end
    
    par_SV = transform_param_ss(par_SV_trans,  [cont.data_on,'_back'], cont);

end