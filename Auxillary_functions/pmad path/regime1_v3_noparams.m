t0 = 0; t_inf = inf; M = 101; 
regime =1;

if inttolspec == 1;
    int_tol;
else
    int_tol = 1e-5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert to an initial condition on u. u0 will be a length n vector,
%interpolated from h0.
idx1 = floor(xl0*length(h0)); idx2 = floor(xu0*length(h0)); %these are the indices corresponding to mensici positions from h0
dx = 1/(length(h0 -1)); 
x = dx*(idx1):dx:dx*(idx2);
v = h0(idx1:idx2);
xq = linspace(x(1), x(end), n); 
v0 = interp1(x,v,xq);
u0 = (xu0 - xl0)*v0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sparsity:
S = spdiags(ones(n+2, 7), -3:3, n+2, n+2); 
S(:, n+1:n+2) = ones(n+2,2);
S(:, n-2:n) = ones(n+2,3); %see droplets for sparsity notes.

%Integrate
dz = 1/n;
h = @(t,v) f_out_1(t,v,n, dz, nu,d,s);
options = odeset('RelTol',int_tol, 'AbsTol',int_tol, 'Events', @(t,v) EventsFcn1(t,v,n), 'Jacobian',@(t,v) PMAD_DfDy(h,{t,v},S));
tspan = [t0 t_inf];
%get solution as a sol structure. 
sol =ode15s(@(t,v) f_out_1(t,v,n, dz, nu,d,s), tspan, [u0 xl0 xu0],options); %get solution as a structure array

%now extract solution as M points:
ts = sol.x;
t = ts;
hs = sol.y';

function f = f_out_1(t,v,n, dz, nu, d,s);


%Original PDE:
%
%h_t = (h^3 h_xxxxx)_x on [xl, xu]
%(h_xxxx = 0 on [0, xl] and [xu ,1])
%
%By changing variables according to: 
%Z = (x-xl)/(xu-xl) and u(Z,t) = h(Z,t)*(xu-xl)
%
%obtain the transformed equation:
%
%u_t + (a(u)u_zzzzz + c(u,z))_z = 0
%where a(u) = -u^3/w^9, c(u,z) = -(z*u*w_t/w + u*xl_t/w), w = xu-xl.
%
%Discretize in space, slve with ODE15 in time.

%First compute the "width" of domain using current last two entries of y:
u = zeros(n+6,1);
u(4:n+3) = v(1:n);
xl = v(n+1); xu = v(n+2); w = xu - xl;

%Boundary conditions:
%
%Top Conditions: u_zz = 0, u_zzz = 0, u_zzzz = -nu*w^6/u:
u_top = (-5*u(n) + 21*u(n+1) - 35*u(n+2) + 35*u(n+3))/16; %from Taylor exp.
uz_top = (-23*u(n)+93*u(n+1)-141*u(n+2)+71*u(n+3))/(24*dz);
%taylor expansion forms of boundary conditions:
A = 0; B = 0; C = -2*nu*w^6 *dz^4/u_top;
u(n+4) = 5*A/264 - B/22 + 13*C/264  +(u(n+1) - 13*u(n+2) +23*u(n+3))/11;
u(n+5) = 2*A/33  - B/22 + 17*C/66   +(u(n+1) - 24*u(n+2) +34*u(n+3))/11;
u(n+6) =19*A/132 - B/22 + 221*C/132 +(-10*u(n+1) - 13*u(n+2) +34*u(n+3))/11;


%Bottom conditions:
%(r2 = A, r3 = B);
% r2 = real((1/xl^3)*(uz_bot*xl/w^2 - 2*u_bot/w +s*xl +2));
% r3 = real((uz_bot/w^2 - 3*r2*xl^2 -s)/(2*xl));
% % r2 = (3*u_bot/w - xl*uz_bot/w^2 -3 - 2*s*xl)/(xl^2);
% % r3 = (uz_bot/(w^2) - 2*r2*xl - s)/(3*xl^2);
% hl_xx = 6*r2*xl + 2*r3; %6Axl +2*B;
% hl_xxx = 6*r2;          %6A

u_bot = (-5*u(7) + 21*u(6) - 35*u(5) + 35*u(4))/16;
uz_bot = -(-23*u(7) + 93*u(6) - 141*u(5) + 71*u(4))/(24*dz); %opposite sign to top
bc1 = 6*(xl*uz_bot/w^2 - 2*u_bot/w +s*xl +2)/(xl^2) + 2*(-2*uz_bot*xl/w^2 +6*u_bot/w -4*s*xl -6)/(xl^2); %hl_xx
bc2 = 6*(xl*uz_bot/w^2 - 2*u_bot/w +s*xl +2)/(xl^3); %hl_xxx

A = 48*dz^2*w^3*bc1;
B = 8*dz^3*w^4*bc2;
C = -2*nu*w^6*dz^4/u_bot;

u(1) = 19*A/132 + B/22 + 221*C/132 +(34*u(4) - 13*u(5) - 10*u(6))/11;
u(2) = 2*A/33   + B/22 + 17*C/66   +(34*u(4) - 24*u(5) + u(6))/11;
u(3) = 5*A/264  + B/22 + 13*C/264  +(23*u(4) - 13*u(5) + u(6))/11;


a = @(u,v) 2*(u.*v).^2 ./(u+v);
A1 = zeros(n+1,1);
A1 = a(u(3:n+3), u(4:n+4)); %midpoint cubed.

zf = 0:dz:1;
xu_t = -(1/4)*(u(n+3) + u(n+4))^2 *(u(n+6) - 5*u(n+5) + 10*u(n+4) - 10*u(n+3) + 5*u(n+2) - u(n+1))/dz^5;
xl_t = -(1/4)*(u(3) + u(4))^2 *(u(6) - 5*u(5) +10*u(4) - 10*u(3) + 5*u(2) - u(1))/dz^5; %NB - no factor of w here as these formulas only used for the flux.
w_t = xu_t - xl_t;


Q =     A1.*((u(6:n+6) - 5*u(5:n+5) + 10*u(4:n+4) - 10*u(3:n+3) + 5*u(2:n+2) - u(1:n+1))/dz^5)...
        +zf(:).*(u(3:n+3) + u(4:n+4))*(1/2)*w_t + (1/2)*(u(3:n+3) + u(4:n+4))*xl_t;
   

f = zeros(n+2,1);
f(1:n) = (1/w^9)*(Q(2:n+1) - Q(1:n))/dz;
f(n+1) = -u_bot^2 *(u(6) - 5*u(5) +10*u(4) - 10*u(3) + 5*u(2) - u(1))/(w^8 * dz^5);
f(n+2) = -u_top^2 *(u(n+6) - 5*u(n+5) +10*u(n+4) - 10*u(n+3) + 5*u(n+2) - u(n+1))/(dz^5 * w^8);
%end
end


%events function to stop it when top meniscus reaches 1 or ends touch
function [position,isterminal,direction] = EventsFcn1(t,v,n)
dz=1/n;
uend = (-5*v(n-3) + 21*v(n-2) - 35*v(n-1) + 35*v(n))/16;
uzend = (-23*v(n-3) + 93*v(n-2) - 141*v(n-1) + 71*v(n))/(24*dz);
w= v(n+2) - v(n+1);
%h(1) = hu'*(1-xu) + hu
endwidth = uend/w + (1-v(n+2))*uzend/w^2;
position = [v(n+2)-(1-1e-4); endwidth; v(n+1)-1e-4]; % Condition 1: meniscus reaches 1. Condition 2: end goes thru zero
isterminal = [1;1;1];  % Halt integration 
direction = [0;0;0];   % The zero can be approached from either direction
end