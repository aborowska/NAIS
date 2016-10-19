clear
% Load data
name = 'Nile.csv';
y = load(name);
n = length(y);
t = 1871:1:1970;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters and initial values
a = zeros(n+1,1);     % filtered state
P = zeros(n+1,1);     % filtered state variance
P(1,1) = 10^7;
sigma2_eps = 15099;
sigma2_eta = 1469.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output of Kalman filter
nu = zeros(n,1);            % prediction errors
F = zeros(n,1);             % prediction variance
K = zeros(n,1);             % = P/F, Kalman gain (regression coefficient of alpha on nu)

for ii = 1:n
   nu(ii,1) = y(ii,1) - a(ii,1);
   F(ii,1) = P(ii,1) + sigma2_eps;
   K(ii,1) = P(ii,1)/F(ii,1);
   a(ii+1,1) = a(ii,1) + K(ii,1)*nu(ii,1);
   P(ii+1,1) = P(ii,1)*(1-K(ii,1)) + sigma2_eta;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State smoothing 
alpha_smooth = zeros(n,1);  % smoothed state
V = zeros(n,1);             % smoothed state variance
r =  zeros(n,1);            % smoothing cumulant
N = zeros(n,1);             % smoothing variance cumulant

L = ones(n,1) - K;
for ii = n:-1:2
    r(ii-1,1) = nu(ii,1)/F(ii,1) + L(ii,1)*r(ii,1);
    alpha_smooth(ii,1) = a(ii,1) + P(ii,1)*r(ii-1,1);
    
    N(ii-1,1) = 1/F(ii,1) + N(ii,1)*L(ii,1)^2;
    V(ii,1) = P(ii,1) - N(ii-1,1)*P(ii,1)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disturbance smoothing (observation disturbance and state disturbance)
u = zeros(n,1);                    % smoothing error
eps_smooth = zeros(n,1);           % smoothed observation disturbance
eta_smooth = zeros(n,1);           % smoothed state disturbance

for ii = n:-1:2
    u(ii,1) = nu(ii,1)/F(ii,1) + r(ii,1)*K(ii,1);
end

eps_smooth = sigma2_eps*u;
eta_smooth = sigma2_eta*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unconditional simulation
eta_un = sqrt(sigma2_eta)*randn(n,1);
alpha_un = zeros(n+1,1);
alpha_un(1,1) = y(1,1);

for ii = 1:n
    alpha_un(ii+1,1) = alpha_un(ii,1) + eta_un(ii,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation smoother 
eps_sim = zeros(n,1);
r_sim = zeros(n,1);
C_sim = zeros(n,1);
W_sim = zeros(n,1);
N_sim = zeros(n,1);

for ii = n:-1:1
    C_sim(ii,1) = sigma2_eps - sigma2_eps^2*(1/F(ii,1) + K(ii,1)^2*N_sim(ii,1));
    W_sim(ii,1) = sigma2_eps*(1/F(ii,1) + K(ii,1)*N_sim(ii,1)*L(ii,1));
    if (ii ~= 1)
        N_sim(ii-1,1) = 1/F(ii,1) + (W_sim(ii,1)^2)/C_sim(ii,1) + L(ii,1)^2*N_sim(ii,1);
    end
    d = sqrt(C_sim(ii,1))*randn;
    
    eps_sim(ii,1) = d + sigma2_eps*(nu(ii,1)/F(ii,1)-K(ii,1)*r_sim(ii,1));
    if (ii ~= 1)
        r_sim(ii-1,1) = nu(ii,1)/F(ii,1) + W_sim(ii,1)*d/C_sim(ii,1) + L(ii,1)*r_sim(ii,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated samples
alpha_sim = zeros(n,1);
eta_sim = zeros(n,1);

for ii = 1:n
    alpha_sim(ii,1) = y(ii,1) - eps_sim(ii,1);
end

for ii = 1:(n-1)
    eta_sim(ii,1) = alpha_sim(ii+1,1) - alpha_sim(ii,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');
 
subplot(2,2,1)
plot(t(1,2:n),alpha_smooth(2:n,1));
hold on
plot(t(1,2:n),alpha_un(2:n,1),'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
hold off
title('smoothed state $$\hat{\alpha}_{t}$$ and sample $$\alpha^{(\cdot)}_{t}$$');
set(gca,'xlim',[1871 1970]);
plotTickLatex2D
 
subplot(2,2,2)
plot(t(1,2:n),alpha_smooth(2:n,1));
hold on
plot(t(1,2:n),alpha_sim(2:n,1),'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
hold off
title('smoothed state $$\hat{\alpha}_{t}$$ and sample $$\tilde{\alpha}_{t}$$');
set(gca,'xlim',[1871 1970]);
plotTickLatex2D

subplot(2,2,3)
plot(t(1,2:n),eps_smooth(2:n,1));
hold on
plot(t(1,2:n),eps_sim(2:n,1),'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
plot(t(1,2:n),zeros(n-1,1));
hold off
title('smoothed observation error $$\hat{\varepsilon}_{t}$$ and sample $$\tilde{\varepsilon}_{t}$$');
set(gca,'xlim',[1871 1970]);
plotTickLatex2D

subplot(2,2,4)
plot(t(1,2:n),eta_smooth(2:n,1));
hold on
plot(t(1,2:n),eta_sim(2:n,1),'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);
plot(t(1,2:n),zeros(n-1,1));
hold off
title('smoothed state error $$\hat{\eta}_{t}$$ and sample $$\tilde{\eta}_{t}$$');
set(gca,'xlim',[1871 1970]);
plotTickLatex2D

suptitle('Simulation')
name = 'Fig4';
saveas(gcf,name,'png'); 
