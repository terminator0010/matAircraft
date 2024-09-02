clc;
clear;
close all;

params.V = [0;0;0];

params.EulerAngles_rad = [1;1;1];
params.Alpha_angle_radps = 0;

params.pqr_radps = [0;0;0];
params.pqr_dot = [0;0;0];

params.H = 100;
params.uvw_mps = [0;0;0];

params.Rotation = [0;0;0];


params.BodyRates_radps = [0;0;0];

params.h0_m = 50;
params.H_ft_z = params.h0_m*(-0.3048);
params.t_k = 50;
params.g_mps = 9.80665;
params.P0 = 101325;

params.e = 2.718281;
params.hTropo = 11000;
params.gConst = 1.4;
params.HR = 287;
params.perfGasEq_ref_kgpm3 = 1.225;

params.GroundForcesAndMoments_N_Nm = [0;0;0];
params.AeroForcesAndMoments_N_Nm = [0;0;0;];
params.ThurstForcesAndMoments_N_Nm = [0;0;0];


params.ExtForces_N = [0;0;0];
params.Ext_Moments_Nm = [0;0;0];
params.LoadFactors_N = [0;0;0];

params.LBE = [0 0 0; 0 0 0; 0 0 0];
params.BodyVelocities = [0;0;0];
params.InertialWind_mps = [0;0;0];

params.TAS_mps = 150;
params.TAS_kt = 150*0.5144444;

params.TC = 0;



[Inertia, y_cg, z_cg, nv, nrho, alfaf_deg, xf_m, zf_m, Tmax, S, c, b, CL0, CL_alpha, CL_elev, CL_AlphaDot, CL_q, CD0, CD_alpha, CD_elev, CY_beta, CY_rud, CY_ail, CY_r, CY_p, Cl_beta, Cl_rud, Cl_ail, Cl_r, Cl_p, Cm0, Cm_alpha, Cm_elev, Cm_AlphaDot, Cm_q, Cn_beta, Cn_rud, Cn_ail, Cn_r, Cn_p] = initACFT_low(params)
params.InertiaTensor_kgm2 = Inertia;
params.Mass_KG = ;
params.nv = nv;
%for turbofan
%nv = 0
%np = 0.7
%for turboprop
params.np = nrho;
params.Tmax = Tmax;
params.y_cg = y_cg;
params.z_cg = z_cg;
params.xf_m = xf_m;
params.zf_m = zf_m;
params.S = S;
params.b = b;
params.c = c;
params.CL0 = CL0;
params.CL_alpha = CL_alpha;
params.CL_elev = CL_elev;
params.CL_AlphaDot = CL_AlphaDot;
params.CL_q = CL_q;
params.CD0 = CD0;
params.CD_alpha = CD_alpha;
params.CD_elev = CD_elev;
params.CY_beta = CY_beta;
params.CY_rud = CY_rud;
params.CY_ail = CY_ail;
params.CY_r = CY_r;
params.CY_p = CY_p;
params.Cl_beta = Cl_beta;
params.Cl_rud = Cl_rud;
params.Cl_ail = Cl_ail;
params.Cl_r = Cl_r;
params.Cl_p = Cl_p;
params.Cm0 = Cm0;
params.Cm_alpha = Cm_alpha;
params.Cm_elev = Cm_elev;
params.Cm_AlphaDot = Cm_AlphaDot;
params.Cm_q = Cm_q;
params.Cn_beta = Cn_beta;
params.Cn_rud = Cn_rud;
params.Cn_ail = Cn_ail;
params.Cn_r = Cn_r;
params.Cn_p = Cn_p;





[t1, spSound_mps, pActual, Mach, cdP, CAS, AlphaVelocities_radps, Alpha_angle_radps, Beta_angle_radps, perfGasEq] = IsaAtmo(params);
params.perfGasEq = perfGasEq;
params.Alpha_angle_radps = Alpha_angle_radps;
params.Beta_angle_radps = Beta_angle_radps;

params.t1_k = t1;
params.CAS_kt = CAS*0.5144444;

[InertiaTensor_kgm2, pqr_dot, pqr_radps, BodyRates_radps, I, BodyVelocities, EulerRates_radps2, EulerAngles_rad, LBE, uvw_dot, uvw_mps, pT] = eQMotion(params);
params.EulerRates_deg = EulerRates_radps2*(pi*180);

[EngineForcesAndMoments, Thrust_N, TAS_Vref] = Propulsion(params);

C_stabAxis = Aerodynamics_Coefficient(params);

AeroForcesMoments = AeroforcesandMoments(params);

C_stabAxis = Aerodynamics_Coefficient(params);

[CxB, CyB, CzB, ClB, CmB, CnB, deltaXcg, Mcg] = ConversionStab2AeroCoeffi(params);



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


