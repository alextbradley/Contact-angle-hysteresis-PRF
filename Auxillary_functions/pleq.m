function pleq(dT, nu, xl)
%Plot the equilibria at these parameters values (note xl must be supplied
%here)
%Use R1_findeq.m to return equilibria
s = 0;
[xu, V, ~,p,~,~] = findeq(dT, nu, s, xl);

%define regions
xxl = linspace(0, xl);
xx  = linspace(xl, xu);
xxu = linspace(xu,1);

%anonymous functions for shape:
%wet region
c = (p*xu^3)/6 + s - (p*xl^3)/6; %linear coefficient
d = (p*xl^4)/24 - (p*xl^3*xu)/6 + (p*xu^4)/8 + s*xu + 1; %constant coeff
hw = @(x) p/24 *(x-xu).^4 + c*(x-xu) + d;

%dry upper
hd = @(x) c*(x-xu) + d;

%dry lower
Q1 = 1/6 * p *(xl - xu); %from matching third derivative
Q2 = 1/2*(1/2*p*(xl-xu)^2 - 6*Q1*xl); %from matching second derivatives
hdl = @(x) 1 + s*x + Q2*x.^2 + Q1*x.^3;
%add drop
borderx = [linspace(-hw(xl),hw(xl)),hw(xx), flip(linspace(-hw(xu),hw(xu))), flip(-hw(xx))];
bordery = [xl*ones(1,100), xx, xu*ones(1,100), flip(xx)];
fill(bordery, borderx, 'b')
%plot
hold on
plot(xxl,hdl(xxl),'k','linewidth', 1.5)
plot(xx,hw(xx), 'k','linewidth', 1.5)
plot(xxu,hd(xxu),'k', 'linewidth', 1.5)
plot(xxl,-hdl(xxl),'k','linewidth', 1.5)
plot(xx, -hw(xx),'k','linewidth', 1.5)
plot(xxu,-hd(xxu),'k','linewidth', 1.5)
xlim([0 1]);
xmax = max([1.2, max(hd(xxu))]);
ylim([-xmax,xmax]*1.1)
end