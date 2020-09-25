function [position,isterminal,direction] = EventsFcn(t,v,n, nu, is_pinned, hysparmax, theta_a, theta_r)
%Events function for the contact angle hysteresis problem.

%Note that we include  small tolerances as a safeguard to ensure that we do
%not spuriously trip events functions at boundaries between integrations.

%define relevant quantities
dz = 1/n;
ut = (-5*v(n-3) + 21*v(n-2) - 35*v(n-1) + 35*v(n))/16;
uzt = (-23*v(n-3) + 93*v(n-2) - 141*v(n-1) + 71*v(n))/(24*dz);
xu = v(n+2);
xl = v(n+1);
w = xu-xl;

%compute pressure and meniscus speeds at x+, x- 
pxl = (7*v(1) - 33*v(2) + 62*v(3) - 58*v(4) + 27*v(5) - 5*v(6))/2 /dz^4 /w^5; %pressure using one sided derivative
uxl = (35*v(1) - 35*v(2) + 21*v(3) - 5*v(4))/16 ; %u at Z = 0
pxu = (7*v(n) - 33*v(n-1) + 62*v(n-2) - 58*v(n-3) + 27*v(n-4) - 5*v(n-5))/2/dz^4 /w^5; %presssure at Z = 1
uxu = (35*v(n) - 35*v(n-1) + 21*v(n-2) - 5*v(n-3));
f = fout(t,v,n, nu, is_pinned, hysparmax); %easier to call fout again, rather than recompute derivatives etc
xl_dot = f(n+1);
xu_dot = f(n+2);
tol = 1e-3;

% Conditions at x = x- 
if is_pinned(1) == 0 %pinned at x-
    %determine contact angle at theta_-
    cosdtheta_l = -pxl*uxl*cosd(theta_a)/nu /w; %avoid taking inverse cosines
    %compute theta_l = acosd(cosdtheta_l);
    
    position(1) = cosdtheta_l - (cosd(theta_a) - tol) ; %if cos(theta_-) passes thru cos(theta_a). Add tolerance to avoid tripping this condition at start of cross over
    direction(1) = 0; 
    
    position(2) = cosdtheta_l - (cosd(theta_r) + tol); %when cos(theta_-) passes thru cos(theta_r) (Note different sense because we're using cos(theta_l) rather than theta_l)
    direction(2) = 0;
    
elseif is_pinned(1) == 1 %receding at x-
    position(1) = xl_dot; %meniscus speed thru 0 (sufficiently negative)
    direction(1) = -1; %only trigger if the speed DECREASES thru 0 (receding associated with positive velocity)
    
    position(2) = 1; %add artificial event function which never goes thru zero so always same number of events conditions
    position(2) = 1;
   
elseif is_pinned(1) == 2 %advancing at x-
    position(1) = xl_dot; %meniscus speed thru 0
    direction(1) = 1; %only trigger if the speed INCREASES thru 0 (advancing associated with negative velocity)
    
    position(2) = 1; %add artificial event function which never goes thru zero so always same number of events conditions
    position(2) = 1;
else
    error('Check pinning conditions: must be a 2 x 1 array with entries in [0,1,2]')
end

%Conditions at x = x+
if is_pinned(2) == 0 %pinned at x+
    %determine contact angle at theta_+
    cosdtheta_u = -pxu*uxu*cosd(theta_a)/nu /w; %avoid taking inverse cosines
    
    position(3) = cosdtheta_u - (cosd(theta_a) - tol) ; %if cos(theta_+) passes thru cos(theta_a). Add tolerance to avoid tripping this condition at start of cross over
    direction(3) = 0; 
    
    position(4) = cosdtheta_u - (cosd(theta_r) + tol); %when cos(theta_+) passes thru cos(theta_r) (Note different sense because we're using cos(theta_l) rather than theta_l)
    direction(4) = 0;
    
    
elseif is_pinned(2) == 1 %receding at x+
    position(3) = xu_dot; %meniscus speed thru 0
    direction(3) = -1; %only trigger if the speed INCREASE thru 0 (receding associated with negative velocity)
    
    position(4) = 1; %add artificial event function which never goes thru zero so always same number of events conditions
    position(4) = 1;
   
elseif is_pinned(2) == 2 %advancing at x+
    position(3) = xu_dot + tol; %meniscus speed thru 0
    direction(3) = -1; %only trigger if the speed DECREASES thru 0 (advancing associated with negative velocity)

    position(4) = 1;
    direction(4) = 1;
else
    error('Check pinning conditions: must be a 2 x 1 array with entries in [0,1,2]')
end

%Conditions independent of the pinning conditions
h1 = uzt*(1-xu)/w^2 + ut/w; %separation at x=1;
position(5) = xu - (1-1e-3); %x+ - 1 increases thru 0
direction(5) = 1;

position(6) = xl - 1e-3; %x- decreases thru 0
direction(6) = -1;  

position(7) = h1; %channel width at x = 1 decreases thru 0  
direction(7) = -1; 

isterminal = ones(length(position), 1); %always terminate the solution
 