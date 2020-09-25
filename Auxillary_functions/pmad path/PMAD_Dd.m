function der = PMAD_Dd(f,args,d,fast,varargin)
%PMAD_Dd Compute directional derivative of a function F(Y) or F(X,Y). 
% DER = PMAD_Dd(F,Y,D) computes the directional derivative of a function
% F(Y). The derivative is returned as a column vector DER. F is a function
% handle. Y and F(Y) are column vectors. Input Y as ARGS. If the function 
% has the form F(X,Y), input the scalar X and vector Y as ARGS = {X,Y}.
%
% If the function is of restricted form and auxiliary functions are 
% used, setting FAST to TRUE can reduce the run time significantly--
% the documentation for details.  The default is FALSE.
%
% DER = PMAD_DFDY(F,ARGS,D,fast,P1,...,PK) is used if the function has the 
% form F(Y,P1,P2,...,PK) or F(X,Y,P1,P2,...PK).  
%
% PMAD_Dd IS LIMITED TO FUNCTIONS THAT DO NOT INVOLVE COMPLEX NUMBERS.

if iscell(args)
    if length(args) == 1
        x = [];
        y = args{1};
    else
        x = args{1};
        y = args{2};
    end
else
    x = [];
    y = args;
end  

nrmd = norm(d);
if nrmd == 0
    error('D must be a non-zero vector.');
else
    d = d(:)/nrmd;
end

if length(y) ~= length(d)
    error('Y and D must have the same length.')
end

if nargin < 4
    fast = false;
end
if nargin < 5
    pars = {};
else
    pars = varargin;
end

% Determine an increment h for complex step differentiation.
h = norm(y); 
if h == 0, h = eps; end
h = eps*h;

Y = y(:);
if ~fast, Y = PMAD_class(Y); end
Y = Y + h*complex(0,1)*d;  

if isempty(x)
    der = imag( feval(f,Y,pars{:}) ) / h;
else
    der = imag( feval(f,x,Y,pars{:}) ) / h;
end

end % PMAD_Dd
      
      
      
      
      