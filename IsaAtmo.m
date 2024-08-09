function [t1, spSound_mps, pActual, Mach, cdP, CAS, AlphaVelocities_radps, Alpha_angle_radps, Beta_angle_radps] = IsaAtmo(params)

  h = params.h0_m;
  t = params.t_k;
  k = -0.0065;
  t0 = 288.15;
  %t_kelvin = t*274,15;
  %h_pes = h*3.2808;
  gamma = params.gConst;
  R = params.HR;
  g = params.g_mps;
  P0 = params.P0;
  e = params.e;
  hTropo = params.hTropo;
  TAS_mps = params.TAS;
  LBE = params.LBE;
  BodyVelocities = params.BodyVelocities;
  InertialWind_mps = params.InertialWind_mps;

  t1 = (t0+k)*h;

  spSound_mps = sqrt(gamma*R*t1);

  pEleven = P0*((t1/t0)^(-g/(k*R)));
  pAbEleven = pEleven^(e^((-g/(R*t1))*(h-hTropo)));

  pActual = pEleven*pAbEleven;

  perfGasEq = pActual*R*t1;


  LEB = transpose(LBE);

  Bodywind_mps = LEB*InertialWind_mps;

  WindSpeed_mps = BodyVelocities+Bodywind_mps;

  TAS_mps = sqrt(WindSpeed_mps(1,1)^2 + WindSpeed_mps(2,1)^2 + WindSpeed_mps(3,1)^2);

  Mach = TAS_mps/spSound_mps;

  cdP = pActual*((((((gamma-1)/2)*Mach^2+1)^(gamma/(gamma-1))))-1);

  Beta_angle_radps = asin(WindSpeed_mps(2,1)/TAS_mps)*(180/pi);

  Alpha_angle_radps = atan2(WindSpeed_mps(3,1),WindSpeed_mps(1,1))*(180/pi);

  AlphaVelocities_radps = diff(Alpha_angle_radps);

  CAS = sqrt((((2*spSound_mps^2)/(gamma-1))*((cdP/P0)+1)^((gamma-1)/gamma))-1);

  end
