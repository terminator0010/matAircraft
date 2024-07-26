function lBE= LBE(params)
  %psi = params.Psi;
  %theta = params.Theta;
  %phi = params.Phi;

  phi = params.EulerAngles_rad(1,1);
  theta = params.EulerAngles_rad(2,1);
  psi = params.EulerAngles_rad(3,1);

  V = params.V;

  %Lpsi = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
  %Ltheta = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
  %Lphi = [1 0 0; 0 cos(phi) sin(ph); 0 -sin(phi) cos(phi)];

  Lpsi = [cos*EulerAngles_rad(psi) sin*EulerAngles_rad(psi) 0; -sin*EulerAngles_rad(psi) cos*EulerAngles_rad(psi) 0; 0 0 1];
  Ltheta = [cos*EulerAngles_rad(theta) 0 -sin*EulerAngles_rad(theta); 0 1 0; sin*EulerAngles_rad(theta) 0 cos*EulerAngles_rad(2,1)];
  Lphi = [1 0 0; 0 cos*EulerAngles_rad(phi) sin*EulerAngles_rad(phi); 0 -sin*EulerAngles_rad(phi) cos*EulerAngles_rad(phi)];

  pV = Lphi*Ltheta*Lpsi;

  LBE = pV*.V;

 % LBE = SICOS_Euler[(5,1)*(6,1) ((1,1)*(2,1)*(6,1))-((4,1)*(3,1)) ((6,1)*(2,1)*(4,1))+((1,1)*u(3,1)); (5,1)*(3,1) (1,1)*(2,1)*(3,1)+(4,1)*(6,1) ((4)*(2)*u(3))-((1)*(6)); -(2); (1,1)*(5,1) (4,1)*(5,1)];


end


