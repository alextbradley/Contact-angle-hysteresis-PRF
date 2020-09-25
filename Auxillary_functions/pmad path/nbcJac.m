function [DgDya,DgDyb,DgDp] = nbcJac(ya,yb,p)
  global bcfun fast
  have_p = (nargin == 3);
  DgDya = PMAD_DfDy(@fya,ya,fast);
  DgDyb = PMAD_DfDy(@fyb,yb,fast);
  if have_p
      DgDp = PMAD_DfDy(@fp,p,fast);
  end

%===Nested functions=========================================
function v = fya(y)
    if have_p
        v = feval(bcfun,y,yb,p);
    else
        v = feval(bcfun,y,yb);
    end
end % fya

function v = fyb(y)
    if have_p
        v = feval(bcfun,ya,y,p);
    else
        v = feval(bcfun,ya,y);
    end
end % fyb

function v = fp(y)
     v = feval(bcfun,ya,yb,y);
end % fp
%==============================================================
  
end % nbcJac  
      
      
      
      