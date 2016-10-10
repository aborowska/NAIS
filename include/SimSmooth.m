% function alpha_cond = SimSmooth(alpha_smooth_obs,par_KFS,RND)
function eps_sim = SimSmooth(S, v, F_inv, eps_smooth, K, L, par_KFS, RND)
%%     % simulation smoothing for a linear Gaussian SSM
% %     % i.e. conditional drawing given the observed sample y_obs
% %     %     c = par_KFS.c; % mean adjustment
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
%  theta=libsimsmo(Z,H,v,FINV,[],K,L,epsilon,RND);
    [n, ~] = size(eps_smooth);
    H = par_KFS.H;
    Z = par_KFS.Z;  
    
    r = zeros(1,S/2);
    eps_bar = zeros(n,S/2);
    eps_sim = zeros(n,S);
    N = 0;
    ind = 1:2:S-1;
    
    for ii = n:-1:1
        C = H(ii,1) - H(ii,1)*(F_inv(ii,1) + N*K(ii,1)^2)*H(ii,1);
        eps_bar(ii,:) = H(ii,1)*(F_inv(ii,1)*v(ii,1) - K(ii,1)*r);
        w = sqrt(C)*RND(ii,:);
        eps_sim(ii,ind) = eps_bar(ii,:) + w; 
        eps_sim(ii,ind+1) = 2*eps_smooth(ii,1) - eps_sim(ii,ind);
        
        W = H(ii,1)*(F_inv(ii,1)*Z(ii,1) - K(ii,1)*N*L(ii,1));
        r = Z(ii,1)*F_inv(ii,1)*v(ii,1) - W*w/C + L(ii,1)*r;
        N = Z(ii,1)*F_inv(ii,1)*Z(ii,1) + W*W/C + N*L(ii,1)^2;
    end
end