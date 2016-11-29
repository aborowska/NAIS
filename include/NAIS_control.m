function cont = NAIS_control(data_on, model)
    % NAIS algorithm specification
    cont.M = 20; % number of the Gauss-Hermite nodes
    cont.tol = 0.0001; % convergence tolerance
    cont.tol_C = 1e-5; % tolerance for the variance 
    [cont.GH.z, cont.GH.h]  = hernodes(cont.M);
    cont.GH.z = -cont.GH.z;
    cont.GH.h = cont.GH.h.*exp(0.5*(cont.GH.z).^2);
    cont.iter_max = 20;
    cont.S = 200;  % number of simulated trajectories in IS estimation
    % then:  RND are normal random numbers, used in SimSmooth, for S/2 simulation paths; 
    % 2 columns for each simulation; S/2 simulations are run one column for eta and one for epsilon; 
    
    % Model specification
    if isempty(strfind(model, 't'))
        cont.err = 'n';
    else
        cont.err = 't';
    end
    
    cont.states = str2num(model(end)); %'sv2', 'svt3'%
    if isempty(cont.states)
        cont.states = 1;      
    end
    
    if data_on
        cont.data_on = 'est'; 
    else
        cont.data_on = 'sim';
    end 
    
    % Interface specification
    cont.print = false; % true for printing iteration info in estimation  
end