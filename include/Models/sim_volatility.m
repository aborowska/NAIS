function [y, theta, alpha] = sim_volatility(par_SV, n, cont)
% Simualtion from a multistate univariate SV model
    m = cont.states; 
    
    c = par_SV(1,1);
    phi = diag(par_SV(:,2:2+m-1));
    s2 = diag(par_SV(:,2+m:2+2*m-1));
    P1 = s2./(1-phi.^2);

    if strcmp(cont.err,'t')
        nu = par_SV(1,end);
    end
%     U = ones(n,m);   

    % the states alpha  
    alpha = zeros(n,m);
    alpha(1,:) = randn(1,m)*sqrt(P1);
    
    eta = randn(n,m)*sqrt(s2);
%     T = phi*U;
    
    for ii = 2:n
%         alpha(ii,:) = T(ii-1)*alpha(ii-1,:) + eta(ii,:);
        alpha(ii,:) = alpha(ii-1,:)*phi + eta(ii,:);
    end
    
    % the signal theta (represents the log-volatility)
%     Z = U;
%     theta = c + Z.*alpha;
    Z = ones(m,1); 
    theta = c + alpha*Z;
    
    % the unobserved volatility sigma2
    sigma2 = exp(theta);

    % the time series of log-returns y
	if strcmp(cont.err,'n')
	    y = sqrt(sigma2).*randn(n,1);
    else % strcmp(cont.err,'t')
		y = sqrt(sigma2.*(nu-2)./nu).*trnd(nu, n, 1);
	end
end
