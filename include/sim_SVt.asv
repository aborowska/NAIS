function y = sim_SVt(par_SV,n)
% simple, one-factor univariate SV model
    c = par_SV(1,1);
    phi = par_SV(1,2);
    s2 = par_SV(1,3);
    P1 = s2/(1-phi^2);
    nu = par_SV(1,4);
    
    U = ones(n,1);   

    % the state alpha  
    alpha = zeros(n,1);
    alpha(1,1) = P1*randn;
    
    eta = sqrt(s2)*randn(n,1);
    T = phi*U;
    
    for ii = 2:n
        alpha(ii,1) = T(ii-1)*alpha(ii-1,1) + eta(ii,1);
    end
    
    % the signal theta (represents the log-volatility)
    Z = U;
    theta = c + Z.*alpha;
    
    % the unobserved volatility sigma2
    sigma2 = exp(theta);

    % the time series of log-returns y
    y = sqrt(sigma2.*(nu-2)./nu).*trnd(nu, n, 1);
end
