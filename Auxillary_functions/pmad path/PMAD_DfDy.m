function Jac = PMAD_DfDy(f,args,fast,vectorized,S,varargin)
%PMAD_DFDY Compute Jacobian dF/dY for a function F(Y) or F(X,Y).
% JAC = PMAD_DFDY(F,ARGS) computes numerically the Jacobian dF/dY of
% a function F(Y). The Jacobian is returned as a full matrix, JAC.
% F is a function handle. Y and F(Y) are column vectors. Input Y 
% as ARGS. If the function has the form F(X,Y), input the scalar X 
% and vector Y as ARGS = {X,Y}.  
%
% For complex functions of many components there are several options
% that can reduce the run time significantly.  If the function is of 
% restricted form and auxiliary functions are used, setting FAST to 
% TRUE avoids the cost of the object-oriented implementation.  See
% the documentation for details about using this option.  The default
% is FALSE.  If VECTORIZED is TRUE, F has been coded so that 
% F([Y1 Y2 ...]) returns the array [F(Y1) F(Y2) ...]. Typically the 
% computation is then much faster. The default is FALSE. A sparse 
% Jacobian is computed when S is a non-empty sparse matrix of 0's and 
% 1's.  A value of 0 for S(i,j) means that component i of the function 
% F does not depend on component j of Y, hence DfDy(i,j) = 0). If the 
% sparsity pattern is favorable, this can reduce the run time dramatically.
%
% JAC = PMAD_DFDY(F,ARGS,...,P1,...,PK) is used if the function has the 
% form F(Y,P1,P2,...,PK) or F(X,Y,P1,P2,...PK).  
%
% PMAD_DFDY IS LIMITED TO FUNCTIONS THAT DO NOT INVOLVE COMPLEX NUMBERS.

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
y = y(:);
ny = length(y);
one2ny = (1:ny)';

if nargin < 3
    fast = false;
end
if nargin < 4 || isempty(vectorized)
    vectorized = false;
end
if nargin < 5 
    S = []; 
end
if nargin < 6
    pars = {};
else
    pars = varargin;
end

if isempty(S)
    g = one2ny;
else
    g = colgroup(S);    
end
ng = max(g);
Y = repmat(y,1,ng);
if ~fast, Y = PMAD_class(Y); end
index = (g-1)*ny + one2ny;
  
% Determine an increment p for complex step differentiation.
p = y;  
%===Deal with degenerate cases=============================
nrmy = norm(y); 
if nrmy == 0, nrmy = eps; end
p( y == 0 ) = eps*nrmy;
%===========================================================
p = eps*p;
Y(index) = Y(index) + complex(0,1)*p;  

if vectorized
    if isempty(x)
        Jac = imag( feval(f,Y,pars{:}) );
    else
        Jac = imag( feval(f,x,Y,pars{:}) );
    end
else
    if isempty(x)
        for col = ng:-1:1
            Jac(:,col) = imag( feval(f,Y(:,col),pars{:}) );
        end
    else
        for col = ng:-1:1
            Jac(:,col) = imag( feval(f,x,Y(:,col),pars{:}) );
        end      
    end
end
nf = length(Jac(:,1));
if isempty(S)
    Jac = Jac ./ repmat(p',nf,1); 
else
    [i j] = find(S);
    Jac = sparse(i,j,Jac((g(j)-1)*nf + i),nf,ny);
    Jac = Jac * sparse(one2ny,one2ny,1 ./ p,ny,ny);
end

end % PMAD_DfDy
      
      
      
      
      