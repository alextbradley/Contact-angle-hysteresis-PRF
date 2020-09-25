function sigma = growthrate_numerics(V, hyspar, nu, N)
%Returns the growth rate of a perturbation to the equilibrium at the point
%of parameter space specified by the inputs. Growth rate obtained by
%solving full equations and fitting to exponential curve at early times.
%Note: requires user input to specify where to perform the fit!

%Inputs:
%
% V      :   Scalar
%            Volume of equilibrium configuration
%
% hyspar :   Scalar
%            Hysteresis parameter (lambda in the paper)
%
% nu     :   Scalar
%            Channel bendability
%
% N      :   Scalar
%            Number of grid points used in finite difference scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Equilibrium Details %%%%%%
% Find the equilibrium by performing a search starting from xl = V and
% increasing until an equilibrium with comparable V found
dT = hyspar + 1;  %findeq accepts parameter called dT, eqs slightly cleaner in dT
xl = 1e-2;
Vfound = 0;
while abs(Vfound - V) > 1e-3
    xl = xl + 1e-4;
    [xu, Vfound_new, hend, p, C, D] = findeq(dT, nu, 0, xl); %compute equilibrium. Third argument specifies clamping angle at x = 0 (equal to 0 in paper)  
    if ~isempty(Vfound_new)
        Vfound = Vfound_new; %only update Vfound if we found a legitimate candidate eq
    end
    
    %if we get beyond xl = 1-V, break and return no NaN
    if xl > 1- V
       warning('Did not find an equilibrium at these parameter values')
       sigma = NaN;
       return 
    end
    
end
xl0 = xl;
xu0 = xu; %initial conditions on full numerics

%equilibrium shape
A = 0; %cubic coefficient
B = 0; %quadratic coeffieicnet
he = @(x) p/24*(x-xu0).^4 + A*(x-xu0).^3 + B*(x-xu0).^2 + C*(x-xu0) + D + 1e-2*sin(x); %sinusoidal perturbation to equilibrium

%%%%%% Numerical Solutions %%%%%%
addpath('Auxillary_functions/pmad path'); %need for auto-diff of Jacobian
dx = 1/N;
max_disp = 1e-4; %how far to run numerics from initial condition
%initial conditions
w0 = xu0 - xl0;
y = dx*((1:N)-0.5)'; %half interior points
x0 = w0*y + xl0; %initial spatial grid points
h0 = he(x0);
u0 = h0*w0; %scaled initial condition
is_pinned = [0,2]; %indicates pinned x_- and advancing x_+

%Sparsity Pattern
S = spdiags(ones(N+2, 7), -3:3, N+2, N+2); %each grid point interacts only with itself and three either side and...
S(:, N+1:N+2) = ones(N+2,2);  %every grid point depends on xu_t and xl_t (these enter into Q)
S(:, N-2:N) = ones(N+2,3);

IC = [u0; xl0; xu0];
h = @(t,v) fout(t,v,N, nu, is_pinned);

options = odeset('RelTol',1e-8, 'AbsTol', 1e-8, ...
    'Events', @(t,v) MyEventFcn(t,v, max_disp, xl0, xu0),... %simple event function that only allows solution to run a certain distance
    'Jacobian',@(t,v) PMAD_DfDy(h,{t,v},S));

sol = ode15s(@(t,v) fout(t,v,N, nu, is_pinned), ...
    [0 1e-1], IC,options); %get solution as a structure array. Run for short amount of time also.

%evalute the solution on a regular grid, and evaluate growth rate on first
%few grid points
t = sol.x;
t = linspace(t(1), t(end), 1e3);
xu = deval(sol, t, N+2);  %get meniscus positions on regular grid
disp = abs(xu - xu(end)); %displacement of meniscus from final position

%perform fit on data up to index specified by idx
idx = 100; %fit on second 50 pts
tt = t(idx:2*idx);
dd = disp(idx:2*idx);
sigma = (log(dd(end)) - log(dd(1)))/(tt(end) - tt(1)); %slope of exponential

% %add dashed line to plot with fit
% exp_fit = dd(1)*exp(sigma*(tt-tt(1)));
% figure(10); clf; hold on
% plot(t,disp, 'mo-');
% set(gca, 'YScale', 'log')
% plot(tt,exp_fit, 'k--');

    function [position,isterminal,direction] = MyEventFcn(t,v, max_disp, xl0, xu0)
        position(1) = abs(v(end-1) - xl0 - max_disp);
        position(2) = abs(v(end) - xu0 -  max_disp);
        %To terminate and direction 
        isterminal = ones(length(position),1); %always terminate
        direction = zeros(length(position),1); %either direction
        
    end
end


    
