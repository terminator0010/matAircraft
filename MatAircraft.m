clc;
clear;
close all;

params.V = [0;0;0];

params.EulerAngles_rad = [1;1;1];

params.pqr_radps = [0;0;0];
params.pqr_dot = [0;0;0];

params.H = 100;
params.H_ft = params.h0*(-0.3048);
params.uvw_mps = [0;0;0];

params.Rotation = [0;0;0];

params.Mass_KG = 10;
params.BodyRates_radps = [0;0;0];

params.h0_m = 50;
params.t_k = 50;
params.g_mps = 9.80665;
params.P0 = 101325;

params.e = 2.718281;
params.hTropo = 11000;
params.gConst = 1.4;
params.HR = 287;

params.GroundForcesAndMoments_N_Nm = [0;0;0];
params.AeroForcesAndMoments_N_Nm = [0;0;0;];
params.ThurstForcesAndMoments_N_Nm = [0;0;0];

params.InertiaTensor_kgm2 = [1 0 1;0 1 0;-1 0 1];


params.ExtForces_N = [0;0;0];
params.Ext_Moments_Nm = [0;0;0];
params.LoadFactors_N = [0;0;0];

params.LBE = [0 0 0; 0 0 0; 0 0 0];
params.BodyVelocities = [0;0;0];
params.InertialWind_mps = [0;0;0];

params.TAS = 150;
params.TAS_kt = 150*0.5144444;




[t1, spSound_mps, pActual, Mach, cdP, CAS, AlphaVelocities_radps, Alpha_angle_radps, Beta_angle_radps] = IsaAtmo(params);
params.t1_k = t1;
params.CAS_kt = CAS*0.5144444;

[InertiaTensor_kgm2, pqr_dot, pqr_radps, BodyRates_radps, I, BodyVelocities, EulerRates_radps2, EulerAngles_rad, LBE, uvw_dot, uvw_mps, pT] = eQMotion(params);

disp('InertiaTensor_kgm2');
disp(InertiaTensor_kgm2);
disp('pqr_dot');
disp(pqr_dot);
disp('pqr_radps');
disp(pqr_radps);
disp('I');
disp(I);
disp('BodyRates_radps');
disp(BodyRates_radps);
disp('BodyVelocities');
disp(BodyVelocities);
disp('EulerRates_radps2');
disp(EulerRates_radps2);
disp('EulerAngles_rad');
disp(EulerAngles_rad);
disp('LBE');
disp(LBE);
disp('uvw_dot');
disp(uvw_dot);
disp('uvw_mps');
disp(uvw_mps);
disp('pT');
disp(pT);

%[EulerAngles_rad, BodyVelocities, Extforces_Nm, Ext_Moments_Nm, LoadFactors_N, LBE, PositionInertial_m, GroundSpeed_mps, Gamma_angle_deg, Track_angle_deg] = eQMotion(params);

%params.LBE = LBE;
%params.BodyVelocities = BodyVelocities;


%[EngineForcesAndMoments, Thrust_N] = Propulsion(params);

%C_stabAxis = Aerodynamics_Coefficient(params);

%[t1, spSound_mps, pActual, Mach, cdP, CAS, AlphaVelocities_radps, Alpha_angle_radps, Beta_angle_radps] = IsaAtmo(params);

%{
disp('Actual Temperature');
disp(t1);
disp('');
disp('Speed of Sound /Mps');
disp(spSound_mps);
disp('');
disp('Actual Pressure');
disp(pActual);
disp('');
disp('Mach Number');
disp(Mach);
disp('');
disp('Dynamic Pressure');
disp(cdP);
disp('');
disp('Calculated Airpressure Speed');
disp(CAS);
%}



%{
[pV, Lpsi, Ltheta, Lphi] = LBE(params);

plV.Lpsi = Lpsi;
plV.Ltheta = Ltheta;
plV.Lphi = Lphi;

pAngR = planeAngRotation(params, plV);
pAngR = planeAngRotation(params, plV);
%}



%pT = planeTranslation(plV, params, m);

%iT = inertiaTensor(params, pAngR);

%disp(pV);

 %% Acceleration for heading error
%{
 figure 1
 plot(plV.Lpsi, plV.Ltheta, plV.Lphi)
 xlabel('Velocities "X", "Y", "Z"', 'fontsize', 14);
 ylabel('Acceleration [g]', 'fontsize', 14);
 title('Plane Velocities','fontsize', 14)
 set(gca, 'fontsize', 14, 'ylim', [0 11], 'xlim', [0 11]);
 set(gcf, 'color', 'w');
 grid on

disp(pAngR);


 figure 2
 plot(pAngR)
 xlabel('Angular Velocities', 'fontsize', 14);
 %ylabel('Acceleration [g]', 'fontsize', 14);
 title('Plane Angular Velocities','fontsize', 14)
 set(gca, 'fontsize', 14, 'ylim', [0 11], 'xlim', [0 11]);
 set(gcf, 'color', 'w');
 grid on
 %}


