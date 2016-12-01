function jaco_inv = jacobian_ss(param_trans, cont)
% Jacobian of the parameter tranformation to get standard errors of the orignal parameters
    
    m = cont.states;
    
    phis = param_trans(:,2:2+m-1);
    sigmas = param_trans(:,2+m:2+2*m-1);
    if strcmp(cont.err,'t')
        nus = param_trans(1,2+2*m:end);
    end
    
    jaco_inv = diag([1, (1+exp(phis)).*(1+exp(-phis)), 1./(exp(sigmas))]);

    if strcmp(cont.err, 't')
        jaco_inv = blkdiag(jaco_inv,diag(1./(exp(nus))));
    end
end
 