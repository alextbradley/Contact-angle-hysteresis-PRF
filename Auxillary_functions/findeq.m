function [xu, V, hend, p, C, D] = findeq(dT, nu, s, xl)
%Compute regime 1 equilbria. This code solves the 10th order polynomial
%derived in the 'Regime1Symbolics.m' code for a given xl, and removes any
%unwanted roots.

%Inputs:
%   dT  : hyseresis parameter (note dT = \lambda + 1  in the paper)
%   nu  : surface tension
%   s   : channel slope (=1 in the paper)
%   xl  : specified lower meniscus position

%Outputs
%   xu  : computed (and allowed) upper meniscus positions
%   V   : associated volumes
%   hend: associated channel ends
%   p0  : associated drop pressures

%polynomial coefficients (from symbolics)
q0 =  nu.*xl.*(s.*xl.^3+xl.^2.*3.0).^2.*-2.4e1-nu.*(dT.*xl.^3+xl.^3.*2.0).*(nu.*xl.^6.*2.0-s.*xl.^3.*2.4e1-xl.^2.*7.2e1+dT.*nu.*xl.^6);
q1 = nu.*(s.*xl.^3+xl.^2.*3.0).^2.*2.4e1-nu.*xl.^2.*(s.*xl.^3+xl.^2.*3.0).*2.88e2+nu.*(dT.*xl.^2.*3.0-xl.^2.*6.0).*(nu.*xl.^6.*2.0-s.*xl.^3.*2.4e1-xl.^2.*7.2e1+dT.*nu.*xl.^6)+nu.*(dT.*xl.^3+xl.^3.*2.0).*(xl.*1.44e2+s.^2.*xl.^3.*2.4e1+s.*xl.^2.*7.2e1+dT.*nu.*xl.^5.*6.0);
q2 = nu.*xl.*((s.*xl.^3+xl.^2.*3.0).*(s.*xl.*3.0+3.0).*2.0+xl.^2.*3.6e1).*-2.4e1+nu.*(dT.*xl.^3+xl.^3.*2.0).*(s.*xl.*2.16e2+nu.*xl.^4.*2.4e1-dT.*nu.*xl.^4.*3.0+7.2e1)+nu.*xl.*(s.*xl.^3+xl.^2.*3.0).*2.88e2-nu.*(dT.*xl.^2.*3.0-xl.^2.*6.0).*(xl.*1.44e2+s.^2.*xl.^3.*2.4e1+s.*xl.^2.*7.2e1+dT.*nu.*xl.^5.*6.0)+dT.*nu.*xl.*(nu.*xl.^6.*2.0-s.*xl.^3.*2.4e1-xl.^2.*7.2e1+dT.*nu.*xl.^6).*3.0;
q3 = nu.*((s.*xl.^3+xl.^2.*3.0).*(s.*xl.*3.0+3.0).*2.0+xl.^2.*3.6e1).*2.4e1+dT.*nu.*(nu.*xl.^6.*2.0-s.*xl.^3.*2.4e1-xl.^2.*7.2e1+dT.*nu.*xl.^6).*3.0-nu.*xl.^2.*(s.*xl.*3.0+3.0).*2.88e2-nu.*(dT.*xl.^2.*3.0-xl.^2.*6.0).*(s.*xl.*2.16e2+nu.*xl.^4.*2.4e1-dT.*nu.*xl.^4.*3.0+7.2e1)+nu.*(dT.*xl.^3+xl.^3.*2.0).*(s.*7.2e1+nu.*xl.^3.*2.4e1+s.^2.*xl.*7.2e1-dT.*nu.*xl.^3.*1.2e1)-dT.*nu.*xl.*(xl.*1.44e2+s.^2.*xl.^3.*2.4e1+s.*xl.^2.*7.2e1+dT.*nu.*xl.^5.*6.0).*3.0;
q4 = nu.*(nu.*xl.^2.*1.8e1-dT.*nu.*xl.^2.*2.7e1).*(dT.*xl.^3+xl.^3.*2.0)-dT.*nu.*(xl.*1.44e2+s.^2.*xl.^3.*2.4e1+s.*xl.^2.*7.2e1+dT.*nu.*xl.^5.*6.0).*3.0-nu.*xl.*(s.*xl.*3.0+3.0).^2.*2.4e1-nu.*(dT.*xl.^2.*3.0-xl.^2.*6.0).*(s.*7.2e1+nu.*xl.^3.*2.4e1+s.^2.*xl.*7.2e1-dT.*nu.*xl.^3.*1.2e1)+nu.*xl.*(s.*xl.*3.0+3.0).*2.88e2-dT.*nu.*xl.*(s.*xl.*2.16e2+nu.*xl.^4.*2.4e1-dT.*nu.*xl.^4.*3.0+7.2e1).*3.0;
q5 = nu.*(s.*xl.*3.0+3.0).^2.*2.4e1-nu.*(nu.*xl.^2.*1.8e1-dT.*nu.*xl.^2.*2.7e1).*(dT.*xl.^2.*3.0-xl.^2.*6.0)-dT.*nu.*(s.*xl.*2.16e2+nu.*xl.^4.*2.4e1-dT.*nu.*xl.^4.*3.0+7.2e1).*3.0-dT.*nu.*xl.*(s.*7.2e1+nu.*xl.^3.*2.4e1+s.^2.*xl.*7.2e1-dT.*nu.*xl.^3.*1.2e1).*3.0-dT.*nu.^2.*xl.*(dT.*xl.^3+xl.^3.*2.0).*1.8e1;
q6 = dT.*nu.*(s.*7.2e1+nu.*xl.^3.*2.4e1+s.^2.*xl.*7.2e1-dT.*nu.*xl.^3.*1.2e1).*-3.0-dT.*nu.^2.*(dT.*xl.^3+xl.^3.*2.0).*9.0+dT.*nu.^2.*xl.*(dT.*xl.^2.*3.0-xl.^2.*6.0).*1.8e1-dT.*nu.*xl.*(nu.*xl.^2.*1.8e1-dT.*nu.*xl.^2.*2.7e1).*3.0;
q7 = dT.*nu.*(nu.*xl.^2.*1.8e1-dT.*nu.*xl.^2.*2.7e1).*-3.0+dT.*nu.^2.*(dT.*xl.^2.*3.0-xl.^2.*6.0).*9.0+dT.^2.*nu.^2.*xl.^2.*5.4e1;
q8 = dT.^2.*nu.^2.*xl.*8.1e1;
q9 = dT.^2.*nu.^2.*2.7e1;
 
%Find roots
    r = roots([q9,q8,q7,q6,q5,q4,q3,q2,q1,q0]);
    
%Remove any unwanted roots
    r  = r(r == real(r)); 
    r = r(r>0);
    r = r(r > xl);
    r = r(r<1);
    
%init sols
    xu = r;
    V = zeros(1,length(xu));
    hend = zeros(1,length(xu));
    p = zeros(1,length(xu));
    C = zeros(1,length(xu));
    D = zeros(1,length(xu));
    
to_store = [];
%Compute associated volumes:
    for i = 1:length(xu)
        %compute pressure
            p(i) = (nu*((dT*(xl^4/24 - (xl^3*xu(i))/6 + xu(i)^4/8))/((xl - xu(i))^4/24 - (xl^3*xu(i))/6 + xl^4/24 + xu(i)^4/8 - (xl - xu(i))*(xl^3/6 - xu(i)^3/6)) - 1))/(s*xu(i) - ((s*xu(i) + s*(xl - xu(i)) + 1)*(xl^4/24 - (xl^3*xu(i))/6 + xu(i)^4/8))/((xl - xu(i))^4/24 - (xl^3*xu(i))/6 + xl^4/24 + xu(i)^4/8 - (xl - xu(i))*(xl^3/6 - xu(i)^3/6)) + 1);
        %and quartic coeffs (again cut and paste)
            c = s-p(i).*xl.^3.*(1.0./6.0)+p(i).*xu(i).^3.*(1.0./6.0);
            d = s.*xu(i)+p(i).*xl.^4.*(1.0./2.4e1)+p(i).*xu(i).^4.*(1.0./8.0)-p(i).*xl.^3.*xu(i).*(1.0./6.0)+1.0;
        %Volume
            V(i) =  -(p(i)/120 *(xl-xu(i)).^5 + c*(xl-xu(i))^2 /2 + d*(xl-xu(i)));
            
        %coefficients
            C(i) = c;
            D(i) = d;
            
        %Channel thickness at menisci
            hu = d;
            hl = p(i)/24 *(xl - xu(i))^4 + c*(xl-xu(i)) + d;
        %Compute errors
            err = max(abs([p(i) + nu/hu, p(i) + nu*dT/hl])); %equations we solved for the polynomial
            
        %Compute free end gap
            h_end = c*(1-xu(i)) + d;
            hend(i) = h_end;
        %Check for solutions (error small, not touching and volume
        %positive)
            if err < 1e-10 && h_end > 0  && V(i) > 0 %&& V(i) < 1
                to_store = [to_store, i];
            end
                
    end
    
 %Remove non-solutions and non-regime 1
    xu = xu(to_store);
    V  = V(to_store);
    hend  = hend(to_store);
    p  = p(to_store);
    C = C(to_store);
    D = D(to_store);
    
    
    if isempty(to_store)
        %warning('no equilibria found')
    end
 end