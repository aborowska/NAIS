% lnP = zeros(n,S);
% for ii = 1:S
%     lnP(:,ii) = log(normpdf(y, 0, sqrt(exp(par_SV.c + theta_sim(:,ii)))));
%      -0.5*(log(2*pi) + log(exp(par_SV.c + theta_sim(:,ii))) - repmat(y.^2),1,S)
% end

% data = csvread('GSPC_ret.csv'); % SP500 log returns
% data = 100*data;
% n = length(data);







%% Output of the Kalmna filter and smoother
figure(4);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,2,1)
plot(y)
hold on
plot(theta_smooth,'r')
hold off
title('Observation $$y_{t}$$ and smoothed signal $$\hat{\theta}_t$$');
legend({'$$y_{t}$$','$$\hat{\theta}_t$$'},'Location','SouthEast','interpreter', 'latex');
plotTickLatex2D

subplot(2,2,2)
plot(V);
title('Smoothed signal variance $$V_t$$');
plotTickLatex2D;

subplot(2,2,3)
plot(nu);
title('Prediction errors $$\nu_t$$');
plotTickLatex2D;

subplot(2,2,4)
plot(F);
title('Prediction variance $$F_t$$');
plotTickLatex2D;

suptitle('Output of Kalman filter and smoother')
