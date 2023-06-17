% Second-order Adams-Bashforth
function [u,f] = ab2(u,f,fm1,dt)
  u = u + 0.5*dt*(3*f - fm1);
end