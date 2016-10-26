function par_NAIS = NAIS_param_copula(par_NAIS, y, par_SC, cont)
    % Algorithm 2: Efficient importance parameters using NAIS
    d = size(par_SC,2);
    n = size(y,1);
    Un = ones(n,1);
      
    M = cont.M;         % number of the Gauss-Hermite nodes
    tol = cont.tol;     % convergence tolerance
    tol_C = cont.tol_C;
    z = cont.GH.z;      % the Gauss-Hermite abscissae z 
    h = cont.GH.h;      % and the associated weights h(z)
    w = sqrt(h);
    iter_max = cont.iter_max;                    

    Um = ones(M,1);
    
    err = 1;
    iter = 0;
    
    if (d == 4)
        nu = par_SC(1,4);
 %         pdf_const = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(nu-2);
    else
        nu = [];
    end
            
    while ((err > tol) && (iter < iter_max))        
        par_KFS = IS_Model(par_NAIS, par_SC);
        b = par_NAIS.b;
        C = par_NAIS.C;
        y_star = b./C;

        [theta_smooth, V_smooth] = KFS(y_star, par_KFS); 
        % compute the smoothed mean theta_smooth
        % and smoothed variance V_smooth based on b, C from previous
        % iteration and the linear SSM using KFS

%         b_new = 0*Un;
%         C_new = 1*Un;
%         for ii = 1:n
% %             % generate the nodes of GH integration
%             theta_GH = theta_smooth(ii,1) + sqrt(V_smooth(ii,1))*z;
%             if cont.link
%                 transf = @(aa) (1 - exp(-aa))./(1 + exp(-aa));
%             else
%                 transf = @(aa) (exp(2*aa)-1)./(exp(2*aa)+1);        
%             end 
%             rho_GH = transf(theta_GH);
% %             % weighted least squares regression
%             if (d == 3)  % if cont.err == 'n'
% % %                 Y = -0.5*(theta_GH + (y(ii,1)^2)./exp(theta_GH)); 
% %                  Y = -0.5*(log(2*pi) + theta_GH + (y(ii,1)^2)./exp(theta_GH)); 
%                 Y = - 0.5*log(1-rho_GH.^2) - 0.5*(y(ii,1).^2 + y(ii,2).^2 - 2*rho_GH.*y(ii,1).*y(ii,2))./(1-rho_GH.^2) ...
%                     + 0.5*y(ii,1).^2 + 0.5*y(ii,2).^2;               
%             else % if cont.err == 't'
%                 Y = log(gamma((nu+2)/2)) + log(gamma(nu/2)) - 2*log(gamma((nu+1)/2)) - 0.5*log(1-rho_GH.^2) ...
%                     - 0.5*(nu+2).*log(1 + (y(ii,1).^2 + y(ii,2).^2 - 2*rho_GH.*y(ii,1).*y(ii,2))./(nu.*(1 - rho_GH.^2))) ...
%                     + 0.5*(nu+1).*log(1 + (y(ii,1).^2)./nu) + 0.5*(nu+1).*log(1 + (y(ii,2).^2)./nu);
%           
% %                  Y = pdf_const - 0.5*(theta_GH + (nu+1)*log(1 + (y(ii,1)^2)./((nu-2).*exp(theta_GH)))); 
% % %                 Y =  - 0.5*(theta_GH + (nu+1)*log(1 + (y(ii,1)^2)./((nu-2).*exp(theta_GH)))); 
%             end
% %             X = [Um, theta_GH, -0.5*theta_GH.^2];
% %             X = repmat(w,1,3).*X;
% %             Y = w.*Y;
% %             beta = (X'*X)\(X'*Y);
% %             b_new(ii,1) = beta(2,1);
% %             if (beta(3,1) < tol_C)
% %                 C_new(ii,1) = tol_C;
% %             else 
% %                 C_new(ii,1) = beta(3,1);        
% %             end
%             
%             beta = EIS_reg(Y,theta_GH,w);
%             b_new(ii,1) = beta(1,1);
%             if (beta(2,1) < tol_C)
%                 C_new(ii,1) = tol_C;
%             else 
%                 C_new(ii,1) = beta(2,1);        
%             end
%              
%         end

        [b_new, C_new] = EIS_reg_vec_copula(y, theta_smooth, V_smooth, ...
            z, w, tol_C, nu, cont.link);
        
        err_b = sum((par_NAIS.b - b_new).^2)/n;
        err_C = sum((par_NAIS.C - C_new).^2)/n;
        err = max(err_b, err_C);
        
        par_NAIS.b = b_new;
        par_NAIS.C = C_new;
        iter = iter + 1;
    end
    if cont.print
        fprintf('NAIS_param iter #: %d.\n', iter)
    end
end