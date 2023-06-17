% Computes center of mass
function Xcm = centroid(body,iks)

  X = body.X;
  Ns = length(X);
  ds = 1/Ns;
  
  x = X(:,1); y = X(:,2);
  xs = real(ifft(iks.*fft(x)));
  ys = real(ifft(iks.*fft(y)));
  A = abs(sum(x.*ys)*ds);
  xcm =  0.5/A*sum((x.^2).*ys)*ds;
  ycm = -0.5/A*sum((y.^2).*xs)*ds;

  Xcm = [xcm,ycm];
  
end



