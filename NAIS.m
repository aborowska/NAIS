function [par_SV, hessian, hessian_tr, theta_smooth] = NAIS(par_SV_init, y, S, cont, RND, options)
    n = length(y);
 
    %%  NAIS algorithm for selecting the importance parameter set [b,C]
    par_init.b = zeros(n,1);           
    par_init.C = ones(n,1);           % V = (Z*P1*Z')';
    par_NAIS = NAIS_param(par_init, y, par_SV_init, cont); % Algorithm 2: Efficient importance parameters via NAIS

    %% MAIN optimisation
    par_SV_init_trans = transform_param_SV(par_SV_init,  [cont.data_on,'_opt']);
%     par_SV_trans = par_SV_init_trans;
    NAIS_max = @(par_SV_trans) NAIS_loglik(par_SV_trans, par_NAIS, y, S, cont, RND);

    if (nargout == 1)
        par_SV_trans = fminunc(NAIS_max, par_SV_init_trans, options);     
    else
        [par_SV_trans,~,~,~,~, hessian]= fminunc(NAIS_max, par_SV_init_trans, options);
%         if strcmp(cont.err,'n')
%             jaco_inv = diag([1, (1-par_SV_trans(1,2))*par_SV_trans(1,2), par_SV_trans(1,3)]);
%         else
%             jaco_inv = diag([1, (1-par_SV_trans(1,2))*par_SV_trans(1,2), par_SV_trans(1,3), par_SV_trans(1,4)-2]);
%         end
%         hessian = jaco_inv*hessian*jaco_inv;
        if strcmp(cont.err,'n')
            jaco_inv = diag([1, (1+exp(par_SV_trans(1,2)))*(1+exp(-par_SV_trans(1,2))), 1/(exp(par_SV_trans(1,3)))]);
        else
            jaco_inv = diag([1, (1+exp(par_SV_trans(1,2)))*(1+exp(-par_SV_trans(1,2))), 1/(exp(par_SV_trans(1,3))), 1/(exp(par_SV_trans(1,4)))]);
        end
        hessian_tr = jaco_inv*(n*hessian)*jaco_inv;
    end    
    
    if (nargout == 4)
        [~, theta_smooth] = NAIS_max(par_SV_trans);
    end
    
    par_SV = transform_param_SV(par_SV_trans,  [cont.data_on,'_back']);

end