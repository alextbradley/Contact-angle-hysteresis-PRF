function w = Cnorm(v,p)
  if nargin == 1
      p = 2;
  end
  if isvector(v)
      switch p
          case 1
              w = sum(Cabs(v));
          case 2
              w = sqrt(sum(Cabs(v).^2));   
          case inf
              w = Cmax(Cabs(v));
          case -inf
              w = Cmin(Cabs(v));
          otherwise
              error('Value of P is not recognized.')
      end
  else
      switch p
          case 1
              w = Cmax(sum(Cabs(v)));
          case 2
              error('In PMAD, NORM(V) and NORM(V,2) are not defined for matrices.')
          case inf
              w = Cmax(sum(Cabs(Ctr(v))));
          case 'fro'
              w = sqrt(sum( sum(Cabs(v).^2) )); 
          otherwise
              error('Value of P is not recognized.')
      end
    
  end
      
 