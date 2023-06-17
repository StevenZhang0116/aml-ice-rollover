% Computes the (scalar) moment of inertia
function I = inertia(body,iks)

  global rho_ice

  X = body.X;
  Xcm = centroid(body,iks);
  Ns = length(X);  
  ds = 1/Ns;

  x = X(:,1)-Xcm(1); y = X(:,2)-Xcm(2);
  ys = real(ifft(iks.*fft(y)));
  I = rho_ice*abs(sum(((1/3)*x.^3 + x.*(y.^2)).*ys))*ds; 

end



