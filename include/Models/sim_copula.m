function [u_sim, theta, rho] = sim_copula(par_true, T, N, link)
    d = size(par_true,2); 

    c = par_true(:,1);
    phi = par_true(:,2);
    sigma2 = par_true(:,3);

    a1 = 0; %c./(1-phi);
    P1 = sigma2./(1-phi.^2);    
    if (d == 4)
        nu = par_true(:,4);
    end
    if link
        transf = @(aa) (1 - exp(-aa))./(1 + exp(-aa));
    else
        transf = @(aa) (exp(2*aa)-1)./(exp(2*aa)+1);        
    end 
    alpha = zeros(N,T);
    u_sim = zeros(T,2,N);   
    
    alpha(:,1) = a1 + sqrt(P1).*randn(N,1);

    for ii = 2:T            
        alpha(:,ii) = phi.*alpha(:,ii-1) + sqrt(sigma2).*randn(N,1);
    end   
    
    theta = c + alpha;
    
    rho = transf(theta);
    if (d == 3)        
        for jj = 1:N
            u_sim(:,:,jj) = bi_copula_n_rnd(rho(jj,:));
        end
    else
        for jj = 1:N
            u_sim(:,:,jj) = bi_copula_t_rnd(rho(jj,:), nu);
        end
    end   

    if (N == 1)
        u_sim = squeeze(u_sim);
    end
    % tic
    % u_sim = zeros(3000,2);
    % for ii=1:T
    %     u_sim(ii,:) = copularnd('t',rho_true(1,ii),mu_true(1,4),1);   
    % end
    % toc%3.895109 

    % tic
    % toc%0.371761
 
end

%  rho_true=[0.99*ones(1,1000),0.5*ones(1,1000),zeros(1,1000)];% 
% figure(1)
% subplot(2,3,1)
% scatter(u_sim(1:1000,1),u_sim(1:1000,2))
% subplot(2,3,2)
% scatter(u_sim(1001:2000,1),u_sim(1001:2000,2))
% subplot(2,3,3)
% scatter(u_sim(2001:3000,1),u_sim(2001:3000,2))
% subplot(2,3,4)
% scatter(u_sim2(1:1000,1),u_sim2(1:1000,2),'r')
% subplot(2,3,5)
% scatter(u_sim2(1001:2000,1),u_sim2(1001:2000,2),'r')
% subplot(2,3,6)
% scatter(u_sim2(2001:3000,1),u_sim2(2001:3000,2),'r')
