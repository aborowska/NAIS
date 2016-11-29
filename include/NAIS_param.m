function par_NAIS = NAIS_param(par_NAIS, y, par_SV, cont)
    % Algorithm 2: Efficient importance parameters using NAIS
    [m, d] = size(par_SV);
    n = length(y);
%     Un = ones(n,1);
      
%     M = cont.M;         % number of the Gauss-Hermite nodes
    tol = cont.tol;     % convergence tolerance
    tol_C = cont.tol_C;
    z = cont.GH.z;      % the Gauss-Hermite abscissae z 
    h = cont.GH.h;      % and the associated weights h(z)
    w = sqrt(h);
    iter_max = cont.iter_max;                    

%     Um = ones(M,1);
    
    err = 1;
    iter = 0;
    
    if strcmp('cont.err','t')
        nu = par_SV(1,4);
        pdf_const = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(nu-2);
    end
            
    while ((err > tol) && (iter < iter_max))        
%         par_KFS = IS_Model(par_NAIS, par_SV);
        par_KFS = IS_Model_multi(par_NAIS, par_SV, cont);
%         b = par_NAIS.b;
%         C = par_NAIS.C;
%         y_star = b./C;
        y_star = par_NAIS.b./par_NAIS.C;
        
%         [theta_smooth, V_smooth] = KFS(y_star, par_KFS);
        [theta_smooth, V_smooth] = KFS_multi(y_star, par_KFS); 

        % compute the smoothed mean theta_smooth
        % and smoothed variance V_smooth based on b, C from previous
        % iteration and the linear SSM using KFS

%% Non-MEXed version
%         b_new = 0*Un;
%         C_new = 1*Un;
%         for ii = 1:n
%             % generate the nodes of GH integration
%             theta_GH = theta_smooth(ii,1) + sqrt(V_smooth(ii,1))*z;
%             % weighted least squares regression
%             if strcmp(cont.err,'n')
% %                 Y = -0.5*(theta_GH + (y(ii,1)^2)./exp(theta_GH)); 
%                 Y = -0.5*(log(2*pi) + theta_GH + (y(ii,1)^2)./exp(theta_GH)); 
%             else % if strcmp(cont.err,'t')
%                 Y = pdf_const - 0.5*(theta_GH + (nu+1)*log(1 + (y(ii,1)^2)./((nu-2).*exp(theta_GH)))); 
% %                 Y =  - 0.5*(theta_GH + (nu+1)*log(1 + (y(ii,1)^2)./((nu-2).*exp(theta_GH)))); 
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

%         [b_new, C_new] = EIS_reg_vec(y, theta_smooth, V_smooth, z, w, tol_C);
%% MEXed vesrion
        if strcmp(cont.err,'n')
            [b_new, C_new] = EIS_reg_vec(y, theta_smooth, V_smooth, z, w, tol_C);
        else % if strcmp(cont.err,'t')
            [b_new, C_new] = EIS_reg_vec_t(y, theta_smooth, V_smooth, z, w, tol_C, nu, pdf_const);
        end
        
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