function pAngR = planeAngRotation(params, plV)

  i = params.I;
  j = params.J;
  k = params.K;

  Lpsi = plV.Lpsi;
  Ltheta = plV.Ltheta;
  Lphi = plV.Lphi;

  pAngR = (Lpsi^-1)*(Ltheta^-1)*(Lphi^-1)*[i;j;k];


end
