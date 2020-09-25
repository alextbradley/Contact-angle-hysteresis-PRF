function [Y,I] = Cmin(X,Z)
% For vectors, MIN(X) is the element in X with the largest real
% part. For matrices, MIN(X) is a row vector containing the element 
% from each column with the largest real part.
%
% [Y,I] = MIN(X) returns the indices of the minimum values in vector I.
% If the values along the first non-singleton dimension contain more
% than one minimal element, the index of the first one is returned.
% 
% MIN(X,Z) returns an array the same size as X and Z with the elements 
% of largest real part taken from X or Z.  Either one can be a scalar.

  if nargin == 1
      if isvector(X)
          [D,I] = min(real(X));
          Y = X(I(1));
      else
          col = size(X,2);
          Y = zeros(1,col);
          [D,I] = min(real(X));
          for m = 1:col
              Y(m) = X(I(m),m);
          end
      end
  else
      sw = real(X) < real(Z);
      Y = sw.*X + (1 - sw).*Z;
  end   
  
      