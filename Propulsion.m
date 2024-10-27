function [EngineForcesAndMoments, Thrust_N, TAS_Vref] = Propulsion(params)

%Generic thrust model
%Thrust = Tc*Tmax*(V/Vref)^nv * (PActual/Pref)^np
%nv = 0 for turbofan and nv = -1 for turboprop
%np = 0,7 - 1.0
%Tc = Throttle Control

%ThrustAngle = AlphaAngle-ThrustForce

%Input
TAS = params.TAS_mps;
perfGasEq = params.perfGasEq;
perfGasEq_ref_kgpm3 = params.perfGasEq_ref_kgpm3;
V_reference_mps = params.V;
nv = params.nv;
np = params.np;
Tc = params.TC;
Tmax = params.Tmax;
xf_m = params.xf_m;
zf_m = params.zf_m;
Alpha_radps = params.Alpha_angle_radps;



TAS_Vref = transpose(TAS/V_reference_mps);
TAS_Vref_nv = TAS_Vref(1,1).^nv;


perfGasEq_Rhoref = (perfGasEq/perfGasEq_ref_kgpm3);
perfGasEq_Rhoref_nrho = perfGasEq_Rhoref(1,1)^perfGasEq;

Thrust_N = (TAS_Vref_nv* perfGasEq_Rhoref_nrho* Tc* Tmax);

Tas_Perf_Gas_Tc_Tmax = (TAS_Vref_nv* perfGasEq_Rhoref_nrho* Tc* Tmax);
x_engine_N = cos(Alpha_radps)* Tas_Perf_Gas_Tc_Tmax;
y_engine_N = 0;
z_engine_N = sin(Alpha_radps)* Tas_Perf_Gas_Tc_Tmax;
Thrust = [x_engine_N; y_engine_N; z_engine_N];
L_engine_m = 0;

xf = Thrust(3,1)*xf_m;
yf = Thrust(2,1)*0;
zf = Thrust(1,1)*zf_m;

xyz_engine = [xf; yf; zf];

pitch_moment_Nm = xyz_engine(1,1) - xyz_engine(3,1);
roll_moment_Nm = 0;
yaw_moment_Nm = 0;

EngineForcesAndMoments = [roll_moment_Nm; pitch_moment_Nm; yaw_moment_Nm];

end
