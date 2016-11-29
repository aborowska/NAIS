% function alpha_cond = SimSmooth(alpha_smooth_obs,par_KFS,RND)
function eps_sim = SimSmooth_multi(S, nu, F_inv, eps_smooth, K, L, par_KFS, RND)
%%     % simulation smoothing for a linear Gaussian SSM
       % de Jong-Shephard method for simulation of disturbances
       % i.e. conditional drawing given the observed sample y_obs
%%     c = par_KFS.c; % mean adjustment
% %         H = par_KFS.H;
% %         Q = par_KFS.Q;
% %         T = par_KFS.T;
% %         R = par_KFS.R;
% %         Z = par_KFS.Z;
% % 
% %         n = length(alpha_smooth_obs);
% % 
% %         % Unconditional simulation
% %         eta_un = sqrt(Q).*RND(:,1);
% %         eps_un = sqrt(H).*RND(:,2);
% %         alpha_un = zeros(n+1,1);
% %         alpha_un(1,1) = alpha_smooth_obs(1,1);
% %         y_un = zeros(n,1);
% %         for ii = 1:n
% %             y_un(ii,1) = Z(ii,1)*alpha_un(ii,1) + eps_un(ii,1);
% %             alpha_un(ii+1,1) = T(ii,1)*alpha_un(ii,1) + R(ii,1)*eta_un(ii,1);
% %         end
% %         alpha_un = alpha_un(1:n,:);
% % 
% %         % Smoothed estimates
% %        alpha_smooth_un = KFS(y_un,par_KFS);
% % 
% %         % Conditional simulation
% %         alpha_cond = alpha_un + alpha_smooth_obs - alpha_smooth_un;

%% conditional simulation of disturbances
    [n, m] = size(K);
    H = par_KFS.H;
    Z = par_KFS.Z;  
    
    r = zeros(S/2,m); % dimensions: for consistency with KFS_multi
    eps_bar = zeros(n,S/2);
    N = zeros(m,m);
    eps_sim = zeros(n,S);
    % antithetics
    ind = 1:2:S-1; 
    
    for ii = n:-1:1
        C = H(ii,1) - H(ii,1)*(F_inv(ii,1) + K(ii,:)*N*K(ii,:)')*H(ii,1);
        eps_bar(ii,:) = H(ii,1)*(nu(ii,1)*F_inv(ii,1) - r*K(ii,:)');
        w = sqrt(C)*RND(ii,:);
        eps_sim(ii,ind) = eps_bar(ii,:) + w; 
        eps_sim(ii,ind+1) = 2*eps_smooth(ii,1) - eps_sim(ii,ind);
        
        L_tmp = squeeze(L(ii,:,:));
        W = H(ii,1)*(Z*F_inv(ii,1) - L_tmp*N*K(ii,:)'); % mx1
        r = bsxfun(@plus, - w'*W'/C + r*L_tmp', nu(ii,1)*F_inv(ii,1)*Z');
        N = Z*F_inv(ii,1)*Z' + W*W'/C + L_tmp*N*L_tmp';
    end
         
end