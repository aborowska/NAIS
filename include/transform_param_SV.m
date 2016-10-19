function par_SV_trans = transform_param_SV(par_SV, mode)
    par_SV_trans = par_SV;
    d = size(par_SV,2);
    
%     if strcmp(mode, 'opt') % transformation for optimization --> unbounded
%         par_SV(1,2) = par_SV(1,2)/2 + 0.5;
%         par_SV_trans(1,2) = log(par_SV(1,2)/(1-par_SV(1,2)));
%         par_SV_trans(1,3) = log(par_SV(1,3));
%         if (d == 4)
%             par_SV_trans(1,4) = log(par_SV(1,4) - 2);
%         end
%     else % if mode == 'back' (transform back)
%         par_SV_trans(1,2) = exp(par_SV(1,2))/(1+exp(par_SV(1,2)));
%         par_SV_trans(1,2) = 2*(par_SV_trans(1,2) - 0.5);
%         par_SV_trans(1,3) = exp(par_SV(1,3));
%         if (d == 4)
%             par_SV_trans(1,4) = 2 + exp(par_SV(1,4));
%         end
%     end
    switch mode
        case 'sim_opt' % transformation for optimization --> unbounded
            par_SV(1,2) = par_SV(1,2)/2 + 0.5;
            par_SV_trans(1,2) = log(par_SV(1,2)/(1-par_SV(1,2)));
            par_SV_trans(1,3) = log(par_SV(1,3));
            if (d == 4)
                par_SV_trans(1,4) = log(par_SV(1,4) - 2);
            end
            
        case 'sim_back' % if mode == 'back' (transform back)
            par_SV_trans(1,2) = exp(par_SV(1,2))/(1+exp(par_SV(1,2)));
            par_SV_trans(1,2) = 2*(par_SV_trans(1,2) - 0.5);
            par_SV_trans(1,3) = exp(par_SV(1,3));
            if (d == 4)
                par_SV_trans(1,4) = 2 + exp(par_SV(1,4));
            end
            
        case 'est_opt'
%             par_SV_trans(1,1) = 100*par_SV(1,1);
            par_SV_trans(1,2) = log(par_SV(1,2)/(1-par_SV(1,2)));
            par_SV_trans(1,3) = log(par_SV(1,3));
            if (d == 4)
                par_SV_trans(1,4) = log(par_SV(1,4) - 2);
            end
            
        case 'est_back'
%             par_SV_trans(1,1) = par_SV(1,1)/100;
            par_SV_trans(1,2) = exp(par_SV(1,2))/(1+exp(par_SV(1,2)));
            par_SV_trans(1,2) = min(0.998, par_SV_trans(1,2));
            par_SV_trans(1,3) = exp(par_SV(1,3));
            if (d == 4)
                par_SV_trans(1,4) = 2 + exp(par_SV(1,4));
                par_SV_trans(1,4) = min(200, par_SV_trans(1,4));
            end
    end
end
           
        