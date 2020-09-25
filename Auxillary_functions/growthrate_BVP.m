function sigma = growthrate_BVP(V, hyspar, nu, init_guess)
%Returns the growth rate sigma of parameter space specified by the inputs
%V, hyspar (= lambda in paper), nu. Growth rate obtained by solving the
%BVP. initial guess at the growth rate must be specified by the argument
%init_guess.

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
% init_gues: Scalar
%            Initial guess at the growth rate sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Equilibrium Details %%%%%%
% Find the equilibrium by performing a search starting from xl = V and
% increasing until an equilibrium with comparable V found
dT = hyspar + 1;  %findeq accepts parameter called dT, eqs slightly cleaner in dT
xl = 1e-2;
Vfound = 0;
while abs(Vfound - V) > 1e-3
    xl = xl + 1e-4;
    [xu, Vfound_new, ~, p, C, D] = findeq(dT, nu, 0, xl); %compute equilibrium. Third argument specifies clamping angle at x = 0 (equal to 0 in paper)
    
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
eq_data = [p, 0, 0, C, D, xu, xl, V]; %array storing information about the equilibrium. Second and third entries indicate no cubic or quadratic terms in equilibrium shape

%%%%%% Solve BVP %%%%%%
sigmainit = init_guess; %initial parameter guess. set perturbation to x_- = 1 wlog
guess = @(x) mat4init(x);

xmesh = linspace(xl, xu);
solinit  =  bvpinit(xmesh, guess, sigmainit); %specify mesh and guess
myodes = @(x,y,  sigma) odes(x,y,nu,eq_data, sigma);
mybcs = @(yxl, yxu, sigma) bcs(yxl, yxu, nu, eq_data, sigma);
options = bvpset('RelTol', 1e-7,'AbsTol', 1e-7, 'NMax', 3000); %really high accuracy
sol = bvp4c(myodes, mybcs,  solinit, options);
pars = sol.parameters;
sigma = pars(1);

    function yinit = mat4init(x)
        %initialising initial guess for BVP4c
        yinit = [sin(x)
            cos(x)
            -sin(x)
            -cos(x)
            sin(x)
            cos(x)];
    end
end



