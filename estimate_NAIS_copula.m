function [par_SV, hessian, hessian_tr, theta_smooth] = estimate_NAIS_copula(par_SC_init, y, fn_jacobian, cont, RND, options)
    n = size(y,1);
 
    %%  NAIS algorithm for selecting the importance parameter set [b,C]
    par_NAIS_init.b = zeros(n,1);           
    par_NAIS_init.C = ones(n,1);           % V = (Z*P1*Z')';
%     par_NAIS = par_NAIS_init;
%     par_SC = par_SC_init;
    
    if (size(par_SC_init,2) == 4)     % For Student's t tranformation depends on the parameter 
        y_init = tinv(y,par_SC_init(1,4));
    else
        y = norminv(y);
        y_init = y;
    end
    par_NAIS = NAIS_param_copula(par_NAIS_init, y_init, par_SC_init, cont); % Algorithm 2: Efficient importance parameters via NAIS

    %% MAIN optimisation
    par_SC_init_trans = transform_param_ss(par_SC_init,  [cont.data_on,'_opt'], cont);
    NAIS_max = @(xx) NAIS_loglik_copula(xx, par_NAIS, y, cont.S, cont, RND);
    
    [par_SC_trans,~,~,~,~, hessian] = fminunc(NAIS_max, par_SC_init_trans, options);
 
    jaco_inv = fn_jacobian(par_SC_trans);
    hessian_tr = jaco_inv*(n*hessian)*jaco_inv; 
    
    if (nargout == 4)
        [~, theta_smooth] = NAIS_max(par_SC_trans);
    end
    
    par_SV = transform_param_ss(par_SC_trans,  [cont.data_on,'_back'], cont);

end