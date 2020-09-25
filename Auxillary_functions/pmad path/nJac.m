function [dFdy,dFdP] = nJac(xa,ya,P)
  global odefun fast vectorized S
  have_p = (nargin == 3);
  dFdy = PMAD_DfDy(@fy,ya,fast,vectorized,S);
  if have_p
      dFdP = PMAD_DfDy(@fp,P,fast,vectorized,S);
  end

%===Nested functions=========================================
function v = fy(y)
    if have_p
        v = feval(odefun,xa,y,P);
    else
        v = feval(odefun,xa,y);
    end
end % fy

function v = fp(y)
     v = feval(odefun,xa,ya,y);
end % fp
%==============================================================
  
end % nJac  
      
      
      
      