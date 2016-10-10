%% SV
% load('results/SML_gspc_updated.mat')
load('results/SML_gspc_updated_short_sv.mat')
signal = theta_smooth;
resid = y./exp(signal/2);

[~,p,jbstat,critval] = jbtest(resid);
figure(1)
qqplot(resid)
kurtosis(resid)

%% SVt
% load('results/SMLt_gspc_updated.mat')
load('results/SMLt_gspc_updated_short_sv.mat')
signal = theta_smooth;
rho = (par_SV_opt(1,4)-2)/par_SV_opt(1,4);
resid = y./(sqrt(rho)*exp(signal/2));

[~,p,jbstat,critval] = jbtest(resid);
figure(2)
qqplot(resid)
kurtosis(resid)