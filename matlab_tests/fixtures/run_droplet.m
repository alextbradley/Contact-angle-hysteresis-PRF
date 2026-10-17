function [vout, tout, t_events_functions, ispinned, istouch] = run_droplet(nu, hysparmax, theta_r, V, xu0, N, tend, logout)
% Test double ONLY. Does not model physics. Forces event 5 (escape).
assert(nargin == 8, 'Generator must provide all eight run_droplet arguments');
assert(nu > 0 && hysparmax > 0 && V > 0 && xu0 > V && N >= 7 && tend > 0);
assert(logout == 0 || logout == 1);
vout = []; tout = []; t_events_functions = []; ispinned = [];
istouch = 5;
end
