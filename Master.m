%Master code for bendotaxis with contact angle hysteresis. Use this code
%for exploring behaviour. Leave default parameters in to produce example
%solutions with same initial conditions and different values of hysparmax
%(lambda_max in the paper). 

%% Preliminaries
clear

%% Parameters and IC
nu        = 1;    %elastocapillary number
hysparmax = 0.02; %hysteresis (dT = cos(theta_r)/cos(theta_a))
theta_r   = 0;  %receding contact angle
theta_a   = acosd(cosd(theta_r)/(1 + hysparmax));
is_pinned = [2,2]; %indicates both menisci advancing (0: pinned, 1: receding, 2: advancing, with first entry corresponding to - meniscus and second to + meniscus)

%Initial Conditions
xu0 = 0.65;     %initial position + meniscus
xl0 = 0.45;     %initial position - meniscus
V = xu0 - xl0; %dimensionless volume
h0 = @(x) 1 + 0*x; %give straight initial condition on h as anonymous function

%Grid and output details
N = 25; %number of grid points
dz = 1/N;
tend = 1; %final time (don't set to infty because trapped droplet solutions with never end!)
u0 = h0(linspace(xl0+dz/2, xu0 - dz/2, N))*V; %initial condition in discretized u space
reltol = 1e-8;
abstol = 1e-8; %use fine tolerances to avoid problems with errors tripping events functions
%% Sparsity and Jacobian
%The problem is du/dt = J*u where J is the Jacobian. J is determined
%using automatic differentiation (external code supplied in pmad path)
auxpath =  [pwd  '/Auxillary_functions/pmad path'];
addpath(auxpath) %Add auxiallary functions.
auxpath =  [pwd  '/Auxillary_functions/'];
addpath(auxpath) %Add auxiallary functions.
%Add here to run the PMAD_DfDy for Jacobian
%Sparsity Pattern
S = spdiags(ones(N+2, 7), -3:3, N+2, N+2); %each grid point interacts only with itself and three either side and...
S(:, N+1:N+2) = ones(N+2,2);  %every grid point depends on xu_t and xl_t (these enter into Q)
S(:, N-2:N) = ones(N+2,3);
    
%% Integrate
nruns = 1;
istouch = 0; %boolean to trigger termination of entire simulation
tnow = 0;
%Integrate
IC = [u0 xl0 xu0];

tout = [];
vout = [];
tic
cases = ["pinned", "receding", "advancing"]; %for updating text output 
while (~istouch && tnow < tend)  
    
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
    tblock = linspace(t(1), t(end), 100);
    v = zeros(length(tblock), N+2);
    for kk = 1:N+2
        v(:, kk) = deval(sol, tblock, kk);
    end

    tout = [tout, tblock];
    vout = [vout;v];
    tnow = t(end);

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
toc
if istouch == 2
    fprintf('Integration ended with t = t_end \n')
elseif istouch == 5
    fprintf('Integration terminated with x_+ = 1 \n');
elseif istouch == 6
    fprintf('Integration terminated with x_- = 0 \n');
elseif istouch == 7
    fprintf('Integration terminated with h(x = 1) = 0 \n');
end
  
%% Check integration with function performing integration 
logout = 1; %set to 1 for logarithmic spacing time points in each block
%[v, t, t_events_functions, ispinned] = run_droplet(nu,hysparmax, theta_r, V, xu0,N, tend, logout);

%% Plots
figure(1); clf; 
plot(tout, vout(:,end), 'r', tout, vout(:,end-1), 'b')

% % % plot as verification
% hold on; plot(t, v(:,end), 'ko', t, v(:,end-1), 'go')
