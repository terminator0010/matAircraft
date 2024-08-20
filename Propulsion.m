function [EngineForcesAndMoments, Thrust_N, TAS_Vref] = Propulsion(params)

%Generic thrust model
%Thrust = Tc*Tmax*(V/Vref)^nv * (PActual/Pref)^np
%nv = 0 for turbofan and nv = -1 for turboprop
%np = 0,7 - 1.0
%Tc = Throttle Control

%ThrustAngle = AlphaAngle-ThrustForce

%Input
params.TAS;
params.rho;
%V_reference_mps
%Rho_reference_kgpm3


TAS_Vref = (TAS/V_reference_mps);
TAS_Vref_nv = TAS_Vref(1,1)^nv;


rho_Rhoref = (rho/Rho_reference_kgpm3);
rho_Rhoref_nrho = rho_Rhoref(1,1)^nrho;

Thrust_N = (TAS_Vref_nv* rho_Rhoref_nrho* Throttle* Tmax);

alfa_radps = (TAS_Vref_nv* rho_Rhoref_nrho* Throttle* Tmax)*rad/180;
x_engine_N = cos(alfa_rapds);
y_engine_N = 0;
z_engine_N = sin(alfa_radps);
Thrust = [x_engine_N; y_engine_N; z_engine_N];
L_engine_m = 0;

xf = Thrust(2,1)*xf_m;
yf = Thrust(2,1)*0;
zf = Thrust(3,1)*zf_m;

xyz_engine = [xf; yf; zf];

pitch_moment_Nm = xyz_engine(1,1) - xyz_engine(1,3);
roll_moment_Nm = 0;
yaw_moment_Nm = 0;

EngineForcesAndMoments = [roll_moment_Nm; pitch_moment_Nm; yaw_moment_Nm];

end
