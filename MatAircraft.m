clc;
clear;
close all;

params.V = 10;
params.R = 0;
params.H = 100;

params.I = 1;
params.J = 1;
params.K = 1;

params.m = 10;

params.h0 = 0;
params.t = 50;
params.g = 9.80665;
params.P0 = 101325;

params.e = 2.718281;
params.hTropo = 11000;
params.gConst = 1.4;
params.HR = 287;

params.TAS = 150;

pIsa = [iT, iSP, iP, iM, icdP, iCAS];

pIsa = IsaAtmo(params);

disp(pIsa);

%{
[pV, Lpsi, Ltheta, Lphi] = planeVelocities(params);

plV.Lpsi = Lpsi;
plV.Ltheta = Ltheta;
plV.Lphi = Lphi;

pAngR = planeAngRotation(params, plV);
pAngR = planeAngRotation(params, plV);
}%



%pT = planeTranslation(plV, params, m);

%pR = planeRotation(params, pAngR);

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



