% Smoothed step function
function val = sigmoid(z,zL,zR)
  global H ds
  val = (zL - zR)./(1 + exp((z - H)/ds)) + zR;
end