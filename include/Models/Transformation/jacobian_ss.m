function jaco_inv = jacobian_sv(param_trans)
% Jacobian of the parameter tranformation to get standard errors of the orignal parameters
    
    d = size(param_trans,2);
    
    if (d == 4)
        jaco_inv = diag([1, (1+exp(param_trans(1,2)))*(1+exp(-param_trans(1,2))), 1/(exp(param_trans(1,3))), 1/(exp(param_trans(1,4)))]);
    else
        jaco_inv = diag([1, (1+exp(param_trans(1,2)))*(1+exp(-param_trans(1,2))), 1/(exp(param_trans(1,3)))]);
    end
end
 