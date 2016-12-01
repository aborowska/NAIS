function par_trans = transform_param_ss(par, mode, cont)  
    m = cont.states;
    
    phis = par(:,2:2+m-1);
    sigmas = par(:,2+m:2+2*m-1);
    if strcmp(cont.err,'t')
        nus = par(1,2+2*m:end);
    end
    
    switch mode
        case 'sim_opt' % transformation for optimization --> unbounded
%             par(1,2) = par(1,2)/2 + 0.5;
            phis_trans = (phis.^2)./cumprod(phis);
            phis_trans = log(phis_trans./(1-phis_trans));
            
            sigmas_trans = log(sigmas);
            if strcmp(cont.err,'t')
                nus_trans = log(nus - 2);
            end
            
        case 'sim_back' % if mode == 'back' (transform back)
            phis_trans = exp(phis)./(1+exp(phis));
%             par_trans(1,2) = 2*(par_trans(1,2) - 0.5);
            phis_trans = cumprod(phis_trans); % for identification: phi1>...>phim
            
            sigmas_trans = exp(sigmas);
            if strcmp(cont.err,'t')
                nus_trans = 2 + exp(nus);
            end
            
        case 'est_opt'
%             par_trans(1,1) = 100*par(1,1);

            phis_trans = (phis.^2)./cumprod(phis);
            phis_trans = log(phis_trans)/(1-phis_trans);
            
            sigmas_trans = log(sigmas);
            if strcmp(cont.err,'t')
                nus_trans = log(nus - 2);
            end
            
        case 'est_back'
%             par_trans(1,1) = par(1,1)/100;

            phis_trans = exp(phis)./(1+exp(phis));
            phis_trans = cumprod(phis_trans);  % for identification: phi1>...>phim         
            phis_trans = min(0.998, phis_trans);
            
            sigmas_trans = exp(sigmas);
            if strcmp(cont.err,'t')
                nus_trans = 2 + exp(nus);
                nus_trans = min(200, nus_trans);
            end
    end
    
    par_trans = [par(1,1), phis_trans, sigmas_trans];   
    if strcmp(cont.err,'t')
       par_trans = [par_trans; nus_trans]; 
    end 
end   