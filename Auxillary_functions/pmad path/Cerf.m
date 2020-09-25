function v = Cerf(z)
zr = real(z);
zi = imag(z);
re = erf(zr);
im = (2/sqrt(pi))*zi .* exp(-zr .^2);
v = complex(re,im);