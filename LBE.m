function lBE= LBE(params)
  psi = params.Psi;
  theta = params.Theta;
  phi = params.Phi;

  V = params.V;

  Lpsi = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
  Ltheta = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
  Lphi = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];

  pV = Lphi*Ltheta*Lpsi;

  lBE = pV*V;

 end
