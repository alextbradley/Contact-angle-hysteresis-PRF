%This code generates the numerics data (see 'numerics_comparison_data.mat')
%used in Figure 11 in manuscript.

%Preliminaries
addpath('Auxillary_functions')
addpath('Auxillary_functions/pmad path'); %need for auto diff of jacobian

%% Main loop
%Parameters
hysparmax_out = 0.02:0.02:0.1; %values of hysteresis parameter to generate data for
nu_out = 1:0.5:6; %above 6 we run risk of trapping (nu \gtrsim 8)
xu0_escape = zeros(length(hysparmax_out), length(nu_out));

for i = 1:length(hysparmax_out)
    hysparmax = hysparmax_out(i);
    
    for j = 1:length(nu_out)
        nu = nu_out(j);
        
        
        %check whether it escapes when close to the base
        xu0 = V + 0.01;
        [vout, tout, t_events_functions, ispinned, istouch] = run_droplet(nu,hysparmax, theta_r, V, xu0,N, tend);
        
        if istouch == 2 %droplet was trapped, so perform the bisection
            fprintf('Droplet is trapped near base, so performing bisection \n')
            xu0_ub = 1; %initial upper bound
            xu0_lb = xu0; %initial lower bound
            xu0_guess = (xu0_lb + xu0_ub)/2; %initial guess
            err = xu0_ub - xu0_lb; %initial error
            
            while err > 5*1e-3 %accuracy of solutions specified by RHS of inequality
                
                
                %perform the integration
                [~, ~, ~, ~, istouch] = run_droplet(nu,hysparmax, theta_r, V, xu0_guess,N, tend);
                
                %update bounds accordingly
                if istouch == 2 %droplet was trapped: update lower bound
                    xu0_lb = xu0_guess;
                elseif istouch == 5
                    xu0_ub = xu0_guess;
                else
                    error('Integration terminated with a condition that is not t = t(end) or x_+ = 1');
                end
                
                %compute error
                err = xu0_ub - xu0_lb;
                xu0_guess = (xu0_lb + xu0_ub)/2;
            end
            xu0_escape(i,j) = xu0_ub;
            
            %account for always stuck
            if xu0_escape(i,j) > 0.95 %droplets squeeze onto end with xu0 = 0.95
                xu0_escape(i,j) = inf;
            end
            
            if xu0_escape(i,j) < V + 0.03
                xu0_escape(i,j) = -1;
            end
            
            
        elseif istouch == 5 %escaped even at this value, so say always escapes
            fprintf('Droplet always escapes with V = %.2f, lambda_max = %.2f, nu = %.2f', V, hysparmax, nu)
            xu0_escape = -1; %set to -1 to identify
            
        else
            error('Initial guess did not terminate with droplet trapped or droplet reaching free end')
        end
        
    end
end

save('Data/numerics_comparison_data.mat', 'hysparmax_out', 'nu_out', 'V','xu0_escape')


