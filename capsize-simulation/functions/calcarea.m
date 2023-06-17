% Computes area
function A = calcarea(body,iks)
    
  X = body.X;
  Ns = length(X);
  ds = 1/Ns;

  x = X(:,1); y = X(:,2);
  ys = real(ifft(iks.*fft(y)));
  A = abs(sum(x.*ys))*ds; 

end

