function Jac = CSDFD(f,args,fast,vectorized,S,varargin)
%CSDFD computes Jacobian by CSD if possible, otherwise FD.
% This function computes a very accurate Jacobian by CSD, complex step
% differentiation, if PMAD is applicable. If it is not, the function 
% resorts to FD, finite differences, which will generally succeed, but
% the Jacobian will be much less accurate.  See the CSD function 
% PMAD_DfDy and the FD function NUMJAC for documentation.

try
    Jac = PMAD_DfDy(f,args,fast,vectorized,S,varargin{:});
catch
    warning('CSD not applicable. Accuracy compromised by using FD.')
    if iscell(args)
        if length(args) == 1
            t = [];
            y = args{1};
        else
            t = args{1};
            y = args{2};
        end
    else
        t = [];
        y = args;
    end
    if nargin < 3
        vectorized = false;
    end
    if nargin < 4
        S = [];
    end
    if nargin < 5
        pars = [];
    else
        pars = varargin;
    end
    thresh = max(norm(y),eps)*eps*ones(size(y));
    if isempty(t)
        fty = f(y);
        t = 0;   % An artificial value that will not be used.
        if isempty(S)
            Jac = numjac(@af,t,y,fty,thresh,[],vectorized); 
        else
            Jac = numjac(@af,t,y,fty,thresh,[],vectorized,S,[]); 
        end
    else
        fty = f(t,y);
        if isempty(S)
            Jac = numjac(@naf,t,y,fty,thresh,[],vectorized); 
        else
            Jac = numjac(@naf,t,y,fty,thresh,[],vectorized,S,[]); 
        end
    end
end

%===Nested functions=======================================================
function v = af(t,y)
    if isempty(pars)
        v = feval(f,y);
    else
        v = feval(f,y,pars{:});
    end
end % af
function v = naf(t,y)
    if isempty(pars)
        v = feval(f,t,y);
    else
        v = feval(f,t,y,pars{:});
    end
end % naf
%==========================================================================

end % CSDFD