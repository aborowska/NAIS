function par_trans = transform_param_ss(par, mode)
    par_trans = par;
    [m,d] = size(par);
    
    switch mode
        case 'sim_opt' % transformation for optimization --> unbounded
%             par(1,2) = par(1,2)/2 + 0.5;
            par_trans(:,2) = (par_trans(:,2).^2)./cumprod(par_trans(:,2));
            par_trans(:,2) = log(par_trans(:,2)./(1-par_trans(:,2)));
            par_trans(:,3) = log(par(:,3));
            if (d == 4)
                par_trans(1,4) = log(par(1,4) - 2);
            end
            
        case 'sim_back' % if mode == 'back' (transform back)
            par_trans(:,2) = exp(par(:,2))./(1+exp(par(:,2)));
%             par_trans(1,2) = 2*(par_trans(1,2) - 0.5);
            par_trans(:,2) = cumprod(par_trans(:,2)); % for identification: phi1>...>phim
            par_trans(:,3) = exp(par(:,3));
            if (d == 4)
                par_trans(1,4) = 2 + exp(par(1,4));
            end
            
        case 'est_opt'
%             par_trans(1,1) = 100*par(1,1);
            par_trans(:,2) = (par_trans(:,2).^2)./cumprod(par_trans(:,2));
            par_trans(:,2) = log(par_trans(:,2)/(1-par_trans(:,2)));
            par_trans(:,3) = log(par(:,3));
            if (d == 4)
                par_trans(1,4) = log(par(1,4) - 2);
            end
            
        case 'est_back'
%             par_trans(1,1) = par(1,1)/100;
            par_trans(:,2) = exp(par(:,2))./(1+exp(par(:,2)));
            par_trans(:,2) = cumprod(par_trans(:,2));  % for identification: phi1>...>phim         
            par_trans(:,2) = min(0.998, par_trans(:,2));
            par_trans(:,3) = exp(par(:,3));
            if (d == 4)
                par_trans(1,4) = 2 + exp(par(1,4));
                par_trans(1,4) = min(200, par_trans(1,4));
            end
    end
end   