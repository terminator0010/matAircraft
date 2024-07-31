function [t1, spSound_mps, pActual, Mach, cdP, CAS] = IsaAtmo(params)

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
  TAS = params.TAS;

  t1 = (t0+k)*h;

  spSound_mps = sqrt(gamma*R*t1);

  pEleven = P0*((t1/t0)^(-g/(k*R)));
  pAbEleven = pEleven^(e^((-g/(R*t1))*(h-hTropo)));

  pActual = pEleven*pAbEleven;

  perfGasEq = pActual*R*t1;

  cdP = pActual*((((((gamma-1)/2)*Mach^2+1)^(gamma/(gamma-1))))-1);

%Alpha_angle = atan2(w/u)
%Beta_angle = sin^-1(v/Vinf)

%Vinf = sqrt(u²+v²+w²)

%Wind_effect
%Input
%LBE
%BodyVelocities
%InertialWind_mps(x,y,z)

%Output
%BodyWind_mps
%WindSpeed_mps
%TrueAirSpeed_mps
%Beta_angle
%Alpha_angle


LEB = transpose(LBE);

Bodywind_mps = LEB*InertialWind_mps;

WindSpeed_mps = BodyVelocities+Bodywind_mps;

TAS_mps = sqrt(WindSpeed_mps(1,1)^2 + WindSpeed_mps(1,2)^2 + WindSpeed_mps(1,3)^2);

Mach = TAS/spSound_mps;

Beta_angle = asin(WindSpeed_mps(1,2)/TASpeed_mps)*(180/pi);

Alpha_angle = atan2(WindSpeed_mps(1,3)/WindSpeed_mps(1,1))*(180/pi);

CAS = sqrt((((2*spSound_mps^2)/(gamma-1))*((cdP/P0)+1)^((gamma-1)/gamma))-1);

  end
