function [iT, iSP, iP, iM, icdP, iCAS] = IsaAtmo(params)

  h = params.h0;
  t = params.t;
  k = -0.0065;
  t0 = 288.15;
  %t_kelvin = t*274,15;
  %h_pes = h*3.2808;
  gamma = params.gConst;
  R = params.HR;
  g = params.g;
  P0 = params.P0;
  e = params.e;
  hTropo = params.hTropo;
  TAS = params.TAS;

  t1 = (t0+k)*h

  spSound_mps = sqrt(gamma*R*t1);

  pEleven = P0*((t1/t0)^(-g/(k*R)));
  pAbEleven = pEleven^(e^((-g/(R*t1))*(h-hTropo)));

  pActual = pEleven*pAbEleven;

  perfGasEq = pActual*R*t1;

  Mach = TAS/spSound_mps;

  cdP = pActual*((((((gamma-1)/2)*Mach^2+1)^(gamma/(gamma-1))))-1);

  CAS = sqrt((((2*spSound_mps^2)/(gamma-1))*((cdP/P0)+1)^((gamma-1)/gamma))-1);

  end
