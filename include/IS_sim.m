% function theta_sim = IS_sim(S,n,theta_smooth,par_KFS,RND)
function theta_sim = IS_sim(S, y_star, eps_smooth, v, F_inv, K, L, par_KFS, RND)

% sample S/2 trajectories for the signal using e.g. the simulation smoother
    eps_sim = SimSmooth(S, v, F_inv, eps_smooth, K, L, par_KFS, RND);
    theta_sim = repmat(y_star,1,S) - eps_sim;
end