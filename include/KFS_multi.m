function [theta_smooth, V, nu, F_inv, eps_smooth, K, L]  = KFS_multi(y, param)
% Univariate multistate Kalman Filter and Smoother with constatnt parameters T, R, Z, Q
% and time varying H (the variance of the linear Gaussian observation)
    n = size(y,1); 
    m = size(param.P1,1);
    
%% Set parameters and initial values
    a = zeros(n+1,m);       % filtered state
    P = zeros(n+1,m,m);     % filtered state variance
    
    P(1,:,:) = param.P1;
    c = param.c;
    H = param.H;
    Q = param.Q;
    d = param.d;
    T = param.T;
    R = param.R;
    Z = param.Z;
    
%% Output of Kalman filter
%%     a: m,n
%     P: m*m,n,
%     nu: 1,n, 
%     F: 1,n 
%     K: m,n 
%     L: m*m,n 
%     invF: 1,n
%            /*a (mx1)*/    
%             /*P (mxm)*/    
%             /*a (mx1)*/                a[i*m+j]=d[i*m+j]+Ta[j]+Kv[j];

%             /*Ta(mx1)=T(mxm)*a(mx1)*/
%             /*Kv(mx1)=K(mxp)*v(px1)*/
%             /*P (mxm)*/
%             /* TPL'(mxm)=TP(mxm)*L'(mxm)*/
%         /* v (px1)*/
%         /* Za(px1)=Z(pxm)*a(mx1)*/
%         /* v(px1)=y(px1)-Za(px1)*/
%         /* F (pxp)*/
%         /* ZP(pxm)=Z(pxm)*P(mxm)*/
%         /* ZPZ'(pxp)=ZP(pxm)*Z'(mxp)*/
%         /* F(pxp)=ZPZ'(pxp)+H*(pxp)*/      
%         /* TP(mxm)=T(mxm)*P(mxm)*/
%         /* ZPT'(pxm)=Z(pxm)*(TP)'(mxm) */ 
% K = ZPR/F
%         /* L(mxp): T-KZ */
%         /* KZ(mxm)=K(mxp)*Z(pxm)*/
%         /* L(mxm): T(mxm)-KZ(mxm) */
%%
    nu = zeros(n,1);            % prediction error
    F = zeros(n,1);             % prediction variance
    F_inv = zeros(n,1);
    K = zeros(n,m);             % = P/F, Kalman gain (regression coefficient of alpha on nu)
    L = zeros(n,m,m);
    
    for ii = 1:n
        nu(ii,1) = y(ii,1) - c - a(ii,:)*Z;
        P_tmp = squeeze(P(ii,:,:));
        F(ii,1) = Z'*P_tmp*Z + H(ii,1);
        F_inv(ii,1) = 1/F(ii,1);
        K(ii,:) = F_inv(ii,1)*Z'*P_tmp*T;
        L(ii,:,:) = T' - Z*K(ii,:);
        a(ii+1,:) = d + a(ii,:)*T + nu(ii,1)*K(ii,:);
        P(ii+1,:,:) = T'*P_tmp*T + R*Q*R' - K(ii,:)'*Z'*P_tmp*T';%K(ii,:)'*F(ii,1)*K(ii,:);     
    end
%     L = T - K.*Z;
    
%% State smoothing 
%%     alpha_smooth = zeros(n,1);  % smoothed state
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
%%    
%% Distrurbance smoothing   
%%     [epsilon,V]=libdistsmo(Z,H,T||v,FINV,K,L,[]);
%     NK=malloc((m)*sizeof(double)); 
%     ZF=malloc((m)*sizeof(double));   
%     ZFv=malloc((m)*sizeof(double));  
%     Lr=malloc((m)*sizeof(double));   
%     ZFZ=malloc((m*m)*sizeof(double));   
%     NL=malloc((m*m)*sizeof(double));   
%     LNL=malloc((m*m)*sizeof(double));  
%     r=malloc((m)*sizeof(double));  
%     N=malloc((m*m)*sizeof(double));  
         
%         /* D (pxp)*/        
%         /* K'NK(pxp)=K'(pxm)*N(mxm)*K(mxp)*/               
%         /* D(pxp)=F^-1(pxp)+KNK(pxp)*/           
%         /* C(pxp)=H(pxp)-HDH(pxp)*/               
%         /* Epsilon (px1)*/       
%         /* Fv(px1)=F(pxp)*v(px1) */
%         /* K'r(px1)=K'(pxm)*r(mx1)*/         
%         /* u(px1)=Fv-Kr*/         
%         /* Z'F(mxp)=Z'(mxp)*F(pxp) */%         
%         /* r (mx1)*/
%         /* ZFv(mx1)=Z'F(mxp)*v(px1) */%            
%         /* L'r(mx1)=L'(mxm)*r(mx1)*/%                     
%         /* N (mxm)*/%         
%         /* Z'FZ(mxm)=Z'F(mx1)*Z(1xm)*/%           
%         /* L'NL(mxm)=L'(mxm)*N(mxm)*L(mxm)*/   
%         /* N(mxm)= Z'FZ(mxm)+LNL(mxm) */
%%  
    
    eps_smooth = zeros(n,1);       % observation error
    C = zeros(n,1);             % observation error variance
    r = zeros(1,m);             % weighted sum of state innovations nu after t
    % u - smoothing error 
    N = zeros(m,m);              % state smoothing error variance
    
    for ii = n:-1:1
        D = F_inv(ii,1) + K(ii,:)*N*K(ii,:)'; % scalar
        C(ii,1) = H(ii,1) - H(ii,1)*D*H(ii,1); % scalar
        u = nu(ii,1)*F_inv(ii) - r*K(ii,:)';     % scalar 
        eps_smooth(ii,1) = H(ii,1)*u;
        r = u*Z' + r*T; %Z'*nu(ii,1)*F_inv(ii,1) + r*squeeze(L(ii,:,:));
        %N = Z*Z'*F_inv(ii,1) + squeeze(L(ii,:,:))*N*squeeze(L(ii,:,:));
        N = Z*D*Z' + T'*N*T - T'*N*K(ii,:)'*Z' - Z*K(ii,:)*N*T; 
    end
    
%% Recovering the smoothed signal: alpha_smooth = y-epsilon;
    theta_smooth = y - eps_smooth;
    V = C;
    
    if (m == 1)
%         P = squeeze(P);
        L = squeeze(L);
    end
end