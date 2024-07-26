clc;
clear;
close all;

params.V[0,0,0];

params.EulerAngles_rad = [phi; theta; psi];

params.P = 1;
params.Q = 1;
params.R = 1;

params.H = 100;

params.I = 1;
params.J = 1;
params.K = 1;

params.Mass_kg = 10;

params.h0_m = 0;
params.t_k = 50;
params.g_mps = 9.80665;
params.P0 = 101325;

params.e = 2.718281;
params.hTropo = 11000;
params.gConst = 1.4;
params.HR = 287;

params.GroundForcesAndMoments_N_Nm = [1;1;1];
params.AeroForcesAndMoments_N_Nm = [1;1;;];
params.ThurstForcesAndMoments_N_Nm = [1;1;1];

params.InertiaTensor_kgm2 = [1 0 1;0 1 0;-1 0 1];

params.EulerAngles_rad;

params.Extforces_N = [1;1;1];
params.Ext_Moments_Nm = [1;1;1];
params.LoadFactors_N = [1;1;1];

params.TAS = 150;

pIsa = [t1, spSound_mps, pActual, Mach, cdP, CAS];

pIsa = IsaAtmo(params);

EulerAngles_rad = bodyVelocities(params);

disp(pIsa);

%{
[pV, Lpsi, Ltheta, Lphi] = LBE(params);

plV.Lpsi = Lpsi;
plV.Ltheta = Ltheta;
plV.Lphi = Lphi;

pAngR = planeAngRotation(params, plV);
pAngR = planeAngRotation(params, plV);
}%



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
 grid on }%



