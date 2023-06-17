% Computes torque due to buoyancy
function T = torque(body,iks)

  global rho_ice rho_air rho_water
    
  X = body.X;
  Xcm = centroid(body,iks);
  Ns = length(X);
  ds = 1/Ns;
  
  x = X(:,1) - Xcm(1); y = X(:,2);
  ys = real(ifft(iks.*fft(y)));
  
  T = 0.5*ds*sum((rho_ice - sigmoid(y,rho_water,rho_air)).*(x.^2).*ys);

end



