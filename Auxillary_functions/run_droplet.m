function [vout, tout, t_events_functions, ispinned, istouch] = run_droplet(nu,hysparmax, theta_r, V, xu0,N, tend,logout)
%This function performs the simulation of droplet dynamics with contact
%angle hysteresis identically to that in the master code. 

% Inputs:
%
% nu        :   Scalar
%               Channel bendability
%
% hysparmax :   Scalar
%               Contact angle hysteresis parameter (= \lambda_{\max} in
%               manuscript)
%
% theta_r   :   Scalar
%               Receding contact angle (= \theta_{+} is manuscript)
%
% V         :   Scalar
%               Dimensionless droplet volume
%
% xu0       :   Scalar
%               Initial position of + meniscus (= x_+^0 in manuscript)
%
% N         :   Positive Integer
%               Number of gridpoints used in discretization
%
% tend      :   Scalar
%               End point of integration (should droplet not reach end)
%
% logout    :   Boolean
%               Specify whether output log (logout = 1) or linear (logout =
%               0) spaced (in each integration block)
%
%
% Outputs:
% 
% vout      :   nruns x (N+2) array 
%               Solution output. nruns is the total number of integrations
%               performed (not determined a priori but rather  by the
%               evolution of the solution -- each change in pinning
%               conditions requires a new run)
%
% tout      :   1 x (nruns x 100) array
%               Time solution output. Each integration has 100 output time
%               points
%
% t_events  :   nruns x 1 array
%               Time points at which each integration is terminated
%
% ispinned  :   nruns x 2 array 
%               Pinning conditions immediately prior to termination of
%               each integrations
% 
% istouch   :   Integer 
%               Specifies conditions on which final integration terminates
%               as follows: istouch = 2: terminates becuase t = tend, 5:
%               terminates because x_+ = 1, 6: terminates because x_- = 0,
%               7: terminates because h(x = 1) reached 0


% Preliminaries
theta_a   = acosd(cosd(theta_r)/(1 + hysparmax)); %advancing contact angle
nruns = 0;
istouch = 0; %boolean to trigger termination of entire simulation
tnow = 0;
cases = ["pinned", "receding", "advancing"]; %for updating text output 

% Initial conditions
is_pinned = [2,2]; %indicates both menisci advancing (0: pinned, 1: receding, 2: advancing, with first entry corresponding to - meniscus and second to + meniscus)
xl0 = xu0 - V;     %initial condition of - meniscus (= x_-^0 in manuscript)
h0 = @(x) 1 + 0*x; %give straight initial condition on h as anonymous function


%Grid and output details
dz = 1/N;
u0 = h0(linspace(xl0+dz/2, xu0 - dz/2, N))*V; %initial condition in discretized u space
IC = [u0 xl0 xu0];
reltol = 1e-8;
abstol = 1e-8; %use fine tolerances to avoid problems with errors tripping events functions

% Sparsity Pattern
S = spdiags(ones(N+2, 7), -3:3, N+2, N+2); %each grid point interacts only with itself and three either side and...
S(:, N+1:N+2) = ones(N+2,2);  %every grid point depends on xu_t and xl_t (these enter into Q)
S(:, N-2:N) = ones(N+2,3);

%Initialize storage arrays
tout = []; %time output
vout = []; %solution output
t_events_functions = [];  %times at which events functions (or end of integration) are triggered
ispinned = []; %pinning conditions at end of each integration


% Perform loop
fprintf("Running droplet dynamics code... \n");
tic
while (~istouch && tnow < tend)
    nruns = nruns + 1;
    
    h = @(t,v) fout(t,v,N,  nu, is_pinned, hysparmax);
    options = odeset('RelTol',reltol, 'AbsTol', abstol, ...
        'Events', @(t,v) EventsFcn(t,v,N,nu, is_pinned,hysparmax, theta_a, theta_r),...
        'Jacobian',@(t,v) PMAD_DfDy(h,{t,v},S));
    
    %perform the integration
    sol = ode15s(@(t,v) fout(t,v,N, nu, is_pinned, hysparmax), ...
        [tnow tend], IC,options); %get solution as a structure array
    
    %extract the solution and add to the array
    t = sol.x;
    v = sol.y;
    IC = v(:,end); %updated initial conditions
    
    %process solution
    if logout
        tblock = logspace(log10(t(1)+ 1e-6), log10(t(end)*0.99999), 300); %deval gets confused with log(0) and can get errors at end points of array 
    else
        tblock = linspace(t(1), t(end), 100);
    end
    v = zeros(length(tblock), N+2);
    for kk = 1:N+2
        v(:, kk) = deval(sol, tblock, kk);
    end
    
    %update the solution
    tout = [tout, tblock];
    vout = [vout;v];
    tnow = t(end);
    t_events_functions = [t_events_functions; tnow]; %time this integration ended
    ispinned = [ispinned; is_pinned]; %add current pinning conditions
    
    
    %check whether termination has occured
    if ~isempty(sol.ie) %if an event triggered terminations
        if any(sol.ie == [5,6,7])
            istouch = sol.ie;
        end
    else %solution terminated with t = tend
        istouch = 2;
    end
    
    %update pinning conditions
    if ~istouch
        is_pinned = update_pinning_conditions(is_pinned, sol);
        txt = join(["Pinning conditions changed:\n x_- now", cases(is_pinned(1) + 1), "\n x_+ now",cases(is_pinned(2)+1) ...
            "\n t = ", num2str(tnow), ...
            "\n x_- = " num2str(vout(end, end-1)), ...
            "\n x_+ = " num2str(vout(end, end)), "\n"]);
        fprintf(txt)
    end
    
end
clock_time = toc;
fprintf('Integration took %.3f secs \n', clock_time);
if istouch == 2
    fprintf('Integration ended with t = t_end \n')
elseif istouch == 5
    fprintf('Integration terminated with x_+ = 1 \n');
elseif istouch == 6
    fprintf('Integration terminated with x_- = 0 \n');
elseif istouch == 7
    fprintf('Integration terminated with h(x = 1) = 0 \n');
end

end

