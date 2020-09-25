function dydx = odes(x,y, nu, eq_data, sigma)
%Returns the derivative vector y = [f', f'',...f''''''] associated 
%with the problem: dy/dx = odes(x,y).
%We linearise about the equilibrium (given by h_e(x) on xl0 < x < xu0) by
%setting h = he(x) + epsilon e^(sigma*t)f(x), 
%        xu = xu0 + epsilon e^(sigma*t)        
%
%Upon linearising, we find the ODE:
%   sigma*f = 1/(3*abs(nu))*(he^3 f_xxxxx)_x
%which must be solved alongside BC (specified in res.m). Constant sigma and
%xl1 are specified by setting d(sigma)/dx = 0, d(xl1)/dx = 0.

% Inputs:
%
%   x   :   (Mx1) array
%           Defines the spatial mesh on which the problem is solved.
%
%   y   :   (6x1) array 
%           y = [f,f',....,f'''''].

%   nu 	:   Scalar
%           Dimensionless surface tension
%
% eq_data: (1x8) array
%          eq_data = [p0, A, B, C, D, xu0, xl0, V]
%          eq_data describes the equilibrium configuration with pressure
%          p0 and shape:
%          he(x) = p0/24 + A*(x-xu)^3 + B*(x-xu)^2 + C*(x-xu) + D.
%          xl0 and xu0 are the equilibrium positions of the menisci.
% 
% sigma :  Scalar
%           Parameter to be determined as part of the solution. Must be
%           passed as an argument.
%
% Outputs:
%
%   dydx :  (6x1) array
%           Returns the derivative dy/dx.


%Extract data from eq_data:
    p0 = eq_data(1);
    A  = eq_data(2);
    B  = eq_data(3);
    C  = eq_data(4);
    D  = eq_data(5);
    xu0 = eq_data(6);
    he  = p0/24 *(x-xu0).^4 + A*(x-xu0).^3 + B*(x-xu0).^2 + C*(x- xu0) + D;
    dhe = p0/24 *4*(x-xu0).^3 + 3*A*(x-xu0).^2 + 2*B*(x-xu0) + C;
    
%Extract:
    f = y(1);
    df = y(2);
    d2f = y(3);
    d3f = y(4);
    d4f = y(5);
    d5f = y(6);
    
 %Derivatives
    dydx = zeros(6,1);
    dydx(1) = df; %f' = df
    dydx(2) = d2f; %f'' = d2f
    dydx(3) = d3f;
    dydx(4) = d4f;
    dydx(5) = d5f;
    dydx(6) = (3*abs(nu).*sigma.*f - 3*he.^2 .*dhe.*d5f)./(he^3); 

end
