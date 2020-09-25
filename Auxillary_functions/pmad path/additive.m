function v = additive(fcn,z,fast)
  if isa(z,'PMAD_class') || fast
      v = [];
      if ~isempty(real(z))
          v = feval(fcn,real(z));
      end
      if ~isempty(imag(z))
          v = v + complex(0,1)*feval(fcn,imag(z));
      end
      if ~fast, v = PMAD_class(v); end
  else
      v = feval(fcn,z);
  end
end % additive
 