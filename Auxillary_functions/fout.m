function f = fout(t,v,N, nu, is_pinned, hysparmax)
%Returns the derivative du/dt = f for the dicretized problem. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRID
% Discretize in space and solve in time using ODE15s. The grid points are
% located at z = dz/2, 3dz/2,...1 - dz/2. Fluxes are evaluated at z = 0,dz,
% 2*dz,...,1-dz, 1. We include three grid points at either end to impose BC
% (the three BC at each end give a system of 3 simultaneous eqs for the
% values of U at these grid point).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%
%   t       :   Scalar
%               Current time of integration
%
%   n       :   Scalar
%               Number of grid points
%
%   v       :   (n + 2) x 1 array
%               Solution at previous timestep. 
%               v = [U(dz/2), U(3*dz/2),...U(1 - dz/2), xl, xu]
%    
%   nu      :   Scalar
%               Elastocapillary number
%
% ispinned  :   (2x1) Integer Array
%               Describes pinning conditions. 0 indicates pinned, 1 indicates
%               receding, 2 indicates advancing. First entry corresponds to x_-,
%               second entry to x_+
% 
%lambdamax  :   Scalar
%               Hysteresis parameter
%
% OUTPUTS
%
%  f_out:   (n + 2) x 1 array
%           f_out = du/dt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialise and extract from previous timestep
dz = 1/N;
U = zeros(N+6,1);
U(4:N+3) = v(1:N);
xl = v(N+1);
xu = v(N+2);
w = xu - xl;

%Determine U, U_z at z = 0 using one sided differences with using grid point
%at z = dz/2, 3*dz/2, 5*dz/2, 7*dz/2
%Using: http://web.media.mit.edu/~crtaylor/calculator.html with sampled
%points at z = 0.5, 1.5, 2.5, 3.5
    Ub = (-5*U(7) + 21*U(6) - 35*U(5) + 35*U(4))/16; %u at z = 0
    Uzb = (23*U(7) - 93*U(6) + 141*U(5) - 71*U(4))/(24*dz); %du/dz at z = 0;
    
%Express U_zz, U_zzz using one sided derviatives
    A0   = 2/xl^2 * (2*xl*Uzb*w   - 3*Ub*w^2 + 3*w^3); % = U_zz
    B0   = 6/xl^3 * (xl*Uzb*w^2 - 2*Ub*w^3 + 2*w^4); % = U_zzz

% Use these to compute U(2), U(3):
U(2) = (B0*dz^3)/2 + 3*A0*dz^2 + 3*U(4) - 2*U(5);
U(3) = (B0*dz^3)/2 + A0*dz^2 + 2*U(4) - U(5);

%Boundary condition at x_- depends on pinning conditions:
if is_pinned(1) == 0 %pinned conditions, impose no flux
    U(1) =  5*U(2) - 10*U(3) + 10*U(4) - 5*U(5) + U(6);
elseif is_pinned(1) == 1 %receding conditions  p = -nu(lambdamax + 1)/h (i.e.  U_zzzz(z = 0) = -nu(lambdamax + 1)*w^6/U)
    U(1) = -nu *(1 + hysparmax)*w^6/Ub *2*dz^4 + 3*U(2) - 2*U(3) - 2*U(4) + 3*U(5) - U(6);
elseif is_pinned(1) == 2 %advancing conditions:  p = -nu/h (i.e.  U_zzzz(z = 0) = -nu*w^6/U)
    U(1) = -nu *w^6 *2*dz^4 / Ub + 3*U(2) - 2*U(3) - 2*U(4) + 3*U(5) - U(6);
else
    error('Check pinning conditions: must be a 2 x 1 array with entries in [0,1,2]')
end

%Top Conditions (no moment or shear)
Ut = (-5*U(N) + 21*U(N+1) - 35*U(N+2) + 35*U(N+3))/16;
Uzt= (-23*U(N)+93*U(N+1)-141*U(N+2)+71*U(N+3))/(24*dz);
U(N+4) =  - U(N+2) + 2*U(N+3);
U(N+5) = - 2*U(N+2) + 3*U(N+3);

%Boundary condition at x_+ depends on pinning conditions:
if is_pinned(2) == 0 %pinned conditions, impose no flux
    U(N+6) =  5*U(N+5) - 10*U(N+4) + 10*U(N+3) - 5*U(N+2) + U(N+1);
elseif is_pinned(2) == 1 %receding conditions  p = -nu(lambdamax + 1)/h (i.e.  U_zzzz(z = 0) = -nu(lambdamax + 1)*w^6/U)
    U(N+6) = -(nu *(1+ hysparmax) * w^6 *2*dz^4 / Ut) - U(N+1) + 3*U(N+2) - 2*U(N+3) - 2*U(N+4) + 3*U(N+5);
elseif is_pinned(2) == 2 %advancing conditions:  p = -nu/h (i.e.  U_zzzz(z = 0) = -nu*w^6/U)
    U(N+6) = -(nu * w^6 *2*dz^4 / Ut) - U(N+1) + 3*U(N+2) - 2*U(N+3) - 2*U(N+4) + 3*U(N+5);
else
    error('Check pinning conditions: must be a 2 x 1 array with entries in [0,1,2]')
end

%Fluxes are prescribed at grid midpoints
a = @(u,v) 2*(u.*v).^2 ./(u+v);
Ucubed = zeros(N+1,1);          %initialize storage
Ucubed = a(U(3:N+3), U(4:N+4)); %midpoint cubed.
Umidpt = (U(3:N+3) + U(4:N+4))/2; %U at midpoints 

zf = (0:dz:1)'; %grid points on which flux evaluated (column vector)
xu_t = -(Ut)^2 *(U(N+6) - 5*U(N+5) + 10*U(N+4) - 10*U(N+3) + 5*U(N+2) - U(N+1))/dz^5/w^8 /3 /abs(nu); %kinematic conditions
xl_t = -(Ub)^2 *(U(6) - 5*U(5) +10*U(4) - 10*U(3) + 5*U(2) - U(1))/dz^5 / w^8 / 3 / abs(nu);
Q    = -Ucubed.*(U(6:N+6) - 5*U(5:N+5) + 10*U(4:N+4) - 10*U(3:N+3) + 5*U(2:N+2) - U(1:N+1))/dz^5 /3 /abs(nu) /w^9 + ...
      -(1-zf).*Umidpt.*xl_t / w - Umidpt.*zf.*xu_t /w;
  
f = zeros(N+2,1);
f(1:N) = -(Q(2:N+1) - Q(1:N))/dz;
if is_pinned(1) == 0
    f(N+1) = 0;
else
    f(N+1) = xl_t;
end
if is_pinned(2) == 0
    f(N+2) = 0;
else
    f(N+2) = xu_t;
end

end

