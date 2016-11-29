function par_KFS = IS_Model_multi(par_NAIS, par_SV, cont)
% Set the parameters for the artificial linear Gaussian model
% To be used in KFS recursions
    
    C = par_NAIS.C;
    m = cont.states; 
    
    c = par_SV(1,1);
    phi = diag(par_SV(:,2:2+m-1));
    sigma2 = diag(par_SV(:,2+m:2+2*m-1));
    P1 = sigma2./(1-phi.^2);
 
    par_KFS.P1 = P1;
%     par_KFS.c = c*Un;
%     par_KFS.H = C.^(-1);
%     par_KFS.Q = sigma2*Un;
%     par_KFS.d = 0*Un;
%     par_KFS.T = phi*Un;
%     par_KFS.R = Un;
%     par_KFS.Z = Un;
    par_KFS.c = c;
    par_KFS.H = C.^(-1);
    par_KFS.Q = sigma2;
    par_KFS.d = zeros(1,m);
    par_KFS.T = phi;
    par_KFS.R = eye(m,m);
    par_KFS.Z = ones(m,1);
end
