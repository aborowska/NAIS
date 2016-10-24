function par_NAIS = NAIS_param2(cont, par_NAIS, y, par_SV)
% Algorithm 2: Efficient importance parameters using NAIS
%     d = size(param_SV,2);
%% Set parameters
    c = par_SV.c;
%     phi = par_SV.phi;
%     s2 = par_SV.sigma2;
    n = par_SV.n;
    P1 = par_SV.P1;
    Un = ones(n,1);
      
    M = cont.M;         % number of the Gauss-Hermite nodes
    tol = cont.tol;     % convergence tolerance
    fast = cont.fast;   % if true, the fast version for weights is run; 
                        % otherwise the full weights are computed
    z = cont.GH.z;      % the Gauss-Hermite abscissae z 
    h = cont.GH.h;      % and the associated weights h(z)
    w = sqrt(h);
    iter_max = cont.iter_max;                    

%     [z, h] = hernodes(M); % the Gauss-Hermite abscissae z and the associated weights h(z)
    Um = ones(M,1);

%% Loop   
    X = repmat(w,1,3).*[ones(M,1),z,-0.5*z.^2];
    X = inv(X'*X)*X'; 
    X = X(2:3,:);
    
    err = 1;
    iter = 0;
    while ((err > tol) && (iter < iter_max))        
        par_KFS = IS_Model(par_NAIS, par_SV);
        b = par_NAIS.b;     % mu
        C = par_NAIS.C;     % V
        b_new = 0*Un;
        C_new = 0*Un;
        
        theta_GH = repmat(b,1,M) + repmat(sqrt(C),1,M).*repmat(z',n,1); % "simulation"
        
        Y = -0.5*(theta_GH + repmat(y.^2,1,M)./(theta_GH.^2));          % logpdfkernel
        Y = Y.*repmat(w',n,1);
        
        % EIS regression
        % [ytilde, H, C] = NAISreg_fast(Y,X,w,mu,V)
        for ii = 1:n
            beta = X*Y(ii,:)';
            b_new(ii,1) = beta(1,1);    
            C_new(ii,1) = beta(2,1);
        end
        % C
        
        y_star = b./C;

        [theta_smooth, V_smooth] = KFS(y_star, par_KFS); 
        % compute the smoothed mean theta_smooth
        % and smoothed variance V_smooth based on b, C from previous
        % iteration and the linear SSM using KFS

        if ~fast
            a = 0.5*(log(abs(C)) - 2*log(pi) - (b.^2)./C);            
        end
        
        b_new = c*Un;
        C_new = P1*Un;
        for ii = 2:n
            % generate the nodes of GH integration
            theta_GH = theta_smooth(ii,1) + sqrt(V_smooth(ii,1))*z;
            % weighted least squares regression
%             p = normpdf(y(ii,1), 0, sqrt(exp(c + theta_GH)));        
%             Y = log(p);
%             Y = -0.5*(log(2*pi) + log(exp(c + theta_GH)) + (y(ii,1)^2)./(exp(c + theta_GH)));
            Y = -0.5*(theta_GH + (y(ii,1)^2)./exp(theta_GH)); %-(1/2)*(theta+y2./exp(theta))
            X = [Um, theta_GH, -0.5*theta_GH.^2];
%             w = exp((z.^2)/2).*h;
            if ~fast
                logg = a(ii,1) + b(ii,1).*theta_GH - 0.5*(theta_GH.^2).*C(ii,1);
%                 wu = exp(log(p) - logg);
                wu = exp(Y - logg);
                w = w.*sqrt(wu);
            end
%             w = sqrt(w);
            X = repmat(w,1,3).*X;
            Y = w.*Y;
            beta = (X'*X)\(X'*Y);
            b_new(ii,1) = beta(2,1);
            C_new(ii,1) = abs(beta(3,1));
            if (C_new(ii,1) < tol)
                C_new(ii,1) = 1;
            end
        end
        err_b = sum((par_NAIS.b - b_new).^2)/n;
        err_C = sum((par_NAIS.C - C_new).^2)/n;
        err = max(err_b, err_C);
        
        par_NAIS.b = b_new;
        par_NAIS.C = C_new;
        iter = iter + 1;
    end
%     fprintf('No. of NAIS_param iterations: %d.\n', iter)
end