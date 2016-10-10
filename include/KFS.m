function [theta_smooth, V, nu, F_inv, epsilon, K, L]  = KFS(y,param)
% Univariate Kalman Filter and Smoother with time varying parameters
% T, R, Z, Q, H
    n = size(y,1); 
      
%% Set parameters and initial values
    a = zeros(n+1,1);     % filtered state
    P = zeros(n+1,1);     % filtered state variance
    
    P(1,1) = param.P1;
    c = param.c;
    H = param.H;
    Q = param.Q;
    d = param.d;
    T = param.T;
    R = param.R;
    Z = param.Z;
    
%% Output of Kalman filter
    nu = zeros(n,1);            % prediction errors
    F = zeros(n,1);             % prediction variance
    F_inv = zeros(n,1);
    K = zeros(n,1);             % = P/F, Kalman gain (regression coefficient of alpha on nu)

    for ii = 1:n
        nu(ii,1) = y(ii,1) - c(ii,1) - Z(ii,1)*a(ii,1);
        F(ii,1) = Z(ii,1)*P(ii,1)*Z(ii,1) + H(ii,1);
        F_inv(ii,1) = 1/F(ii,1);
        K(ii,1) = T(ii,1)*P(ii,1)*Z(ii,1)/F(ii,1);
        a(ii+1,1) = d(ii,1) + T(ii,1)*a(ii,1) + K(ii,1)*nu(ii,1);
        P(ii+1,1) = T(ii,1)*P(ii,1)*T(ii,1) + R(ii,1)*Q(ii,1)*R(ii,1) - K(ii,1)*F(ii,1)*K(ii,1);     
    end
    L = T - K.*Z;
    
%% State smoothing 
%     alpha_smooth = zeros(n,1);  % smoothed state
%     V = zeros(n,1);             % smoothed state variance
%     r = zeros(n,1);             % smoothing cumulant
%     N = zeros(n,1);             % smoothing variance cumulant
% 
%     for ii = n:-1:2
%         r(ii-1,1) = nu(ii,1)/F(ii,1) + L(ii,1)*r(ii,1);
%         alpha_smooth(ii,1) = a(ii,1) + P(ii,1)*r(ii-1,1);
%         N(ii-1,1) = 1/F(ii,1) + L(ii,1)*N(ii,1)*L(ii,1);
%         V(ii,1) = P(ii,1) - P(ii,1)*N(ii-1,1)*P(ii,1);
%     end   
%     theta_smooth = c + Z.*alpha_smooth;
    
%% Distrurbance smoothing   %         [epsilon,V]=libdistsmo(Z,H,T||v,FINV,K,L,[]);
    epsilon = zeros(n,1);       % smoothed disturbance
    C = zeros(n,1);             % smoothed disturbance variance
    N = 0;
    r = 0;
    
    for ii = n:-1:1
        D = F_inv(ii,1) + K(ii,1)*N*K(ii,1);
        C(ii,1) = H(ii,1) - H(ii,1)*D*H(ii,1);
        u = nu(ii)*F_inv(ii) - K(ii)*r;
        epsilon(ii,1) = H(ii,1)*u;
        r = Z(ii,1)*nu(ii,1)*F_inv(ii,1) + L(ii,1)*r;
        N = Z(ii,1)*Z(ii,1)*F_inv(ii,1) + L(ii,1)*N*L(ii,1);
    end
    
%% Recovering the smoothed signal: alpha_smooth = y-epsilon;
    theta_smooth = y - epsilon;
    V = C;
end