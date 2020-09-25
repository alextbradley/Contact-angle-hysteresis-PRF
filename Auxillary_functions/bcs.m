function res = bcs(yxl, yxu, nu, eq_data, sigma)
%Returns the residuals of the boundary conditions from the linearisation.
%
%
%We linearise about the equilibrium (given by h_e(x) on xl0 < x < xu0) by
%setting h = he(x) + epsilon e^(sigma*t)f(x), 
%        xu = xu0 + epsilon e^(sigma*t)          [Setting xu1 = 1 wlog]
%
%The linearised ODEs are described in odes.m
%
%The full boundary conditions are:
%
% At x = xl:
%  i)   h_xx   = (2/xl^2)*(2*xl*h_x - 3*h + 3)
% ii)   h_xxx  = (6/xl^3)*(xl*h_x - 2*h + 2)
%iii)   h_xxxxx = 0
%
% At x = xu
% iv)   h_xxxx = -nu/h
%  v)   h_xx   = 0
% vi)   h_xxx  = 0
%
% And kinematic condition:
%viii)  d(xu)/dt  = -h^2 /(3*abs(nu)) *h_xxxxx at x = xu
%
% The linearised boundary conditions (by inserting the expansion (*) into
% these) are:
% At x = xl0:
% 
% Effective BC:
%  i)    f_xx - 2/xl0^2 *(2*xl0*f_x - 3*f) = 0
% ii)     f_xxx -  6/(xl0^3) *( xl0*f_x -2*f)
%
% No flux BC:
%iii)   f_xxxxx = 0
%
%
% At x = xu0:
% Pressure BC:
%iv)    f_xxxx - nu/he *(xl1*he_x/he + f/h_e) = 0
%v)      f_xx = 0
%vi)     f_xxx - xu1*nu/he = 0
%   
% Kinematic Condition:
%vii)  sigma*f + he^2/(3*abs(nu)) *f_xxxxx at x = xu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs: 
%   yxl     :   (6x1) array
%               Contains values of y = [f,f',....,f''''']
%               evaluated x = xl0.
%
%   yxu     :   (6x1) array
%               Contains values of y = [f,f',....,f'''''] 
%               evaluated at x = xu0.
%
%   nu      :   Scalar
%               Dimensionless surface tension
%
%   eq_data :   (1x8) array
%               eq_data = [p0, A, B, C, D, xu0, xl0, V]
%               eq_data describes the equilibrium configuration with pressure
%               p0 and shape:
%               he(x) = p0/24 + A*(x-xu)^3 + B*(x-xu)^2 + C*(x-xu) + D.
%               xl0 and xu0 are the equilibrium positions of the menisci
%
% sigma     :   Scalar
%               Growth rate that is to be determined as part of the
%               solution.
%
%Outputs:
%
%   res     :   (7x1) array
%               Returns the residuals from the boundary conditions.
%

%Extract data from eq_data:
p0 = eq_data(1);
A  = eq_data(2);
B  = eq_data(3);
C  = eq_data(4);
D  = eq_data(5);
xu0 = eq_data(6);
xl0 = eq_data(7);

he   = @(x) p0/24 *(x-xu0)^4 + A*(x-xu0)^3 + B*(x-xu0)^2 + C*(x- xu0) + D;
dhe  = @(x) p0/24 *4*(x-xu0)^3 + 3*A*(x-xu0)^2 + 2*B*(x-xu0) + C;
d2he = @(x) p0/24 *12*(x-xu0).^2 + 6*A*(x-xu0) + 2*B;
d3he = @(x) p0*(x-xu0) + 6*A;
d4he = @(x) p0;
    
%Boundary conditions at xl0:
f   = yxl(1); %= f(xl0)
df  = yxl(2); %df/dx at x = xl0
d2f = yxl(3); %d2f/dx2 at x = xl0
d3f = yxl(4); % etc
d4f = yxl(5);
d5f = yxl(6); % etc
    
%Residuals at x = xl
res(1) = (xl0)^2 *d2f -2*(- 3*f + 2*xl0*df); %Effective BC1 from h'' = 2/xl^2 *(2*xl*h' - 3*h + 3)
res(2) = (xl0)^3*d3f -6*(xl0*df - 2*f); %Effective BC2 from h''' = 6/xl^3 *(xl*h' - 2h + 2
res(3) = d5f; %No flux BC at x = xl0
    
%Boundary conditions at xu0:
f   = yxu(1); %= f(xl0)
df  = yxu(2); %df/dx at x = xl0
d2f = yxu(3); %d2f/dx2 at x = xl0
d3f = yxu(4); % etc
d4f = yxu(5);
d5f = yxu(6); % etc
    
%Residuals at x = xu
xu1 = 1; %wlog 
res(4) = d2f + xu1*d3he(xu0); %no moment
res(5) = d3f + d4he(xu0)*xu1; %no shear
res(6) = d4f + d4he(xu0)*(dhe(xu0)*xu1/he(xu0) + f/he(xu0));  %Pressure Bc
res(7) = 3*abs(nu) * sigma*xu1 + (he(xu0)^2) *d5f;%Kinematic condition
end 