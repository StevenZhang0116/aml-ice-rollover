% Computes buoyancy force
function Fb = buoyancy(body,iks)

  global rho_ice rho_air rho_water
    
  X = body.X;
  Ns = length(X);
  ds = 1/Ns;

  x = X(:,1); y = X(:,2);
  ys = real(ifft(iks.*fft(y)));
  Fb = sum(x.*ys.*(rho_ice - sigmoid(y,rho_water,rho_air)))*ds;
end



