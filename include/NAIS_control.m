function cont = NAIS_control(data_on)
    cont.M = 20; % number of the Gauss-Hermite nodes
    cont.tol = 0.0001; % convergence tolerance
    cont.tol_C = 1e-5; % tolerance for the variance 
    [cont.GH.z, cont.GH.h]  = hernodes(cont.M);
    cont.GH.z = -cont.GH.z;
    cont.GH.h = cont.GH.h.*exp(0.5*(cont.GH.z).^2);
    cont.iter_max = 20;
    % cont.err = 'n';
    cont.S = 200;  % number of simulated trajectories in IS estimation
    % then:  RND are normal random numbers, used in SimSmooth, for S/2 simulation paths; 
    % 2 columns for each simulation; S/2 simulations are run one column for eta and one for epsilon; 
    if data_on
        cont.data_on = 'est'; 
    else
        cont.data_on = 'sim';
    end 
    cont.print = false; % true for printing iteration info in estimation  
end