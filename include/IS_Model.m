function par_KFS = IS_Model(par_NAIS, par_SV)
% Set the parameters for the artificial linear Gaussian model
% To be used in KFS recursions
    
    C = par_NAIS.C;
    n = length(C);
    Un = ones(n,1);
    
    c = par_SV(1,1);
    phi = par_SV(1,2);
    sigma2 = par_SV(1,3);
    P1 = sigma2/(1-phi^2);
 
    par_KFS.P1 = P1;
    par_KFS.c = c*Un;
    par_KFS.H = C.^(-1);
    par_KFS.Q = sigma2*Un;
    par_KFS.d = 0*Un;
    par_KFS.T = phi*Un;
    par_KFS.R = Un;
    par_KFS.Z = Un;
end
