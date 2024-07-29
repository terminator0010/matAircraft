function EulerAngles_rad, BodyVelocities, EulerAngles_rad, Extforces_N, Ext_Moments_Nm, LoadFactors_N = eQMotion(params)

%Input
ExtForces_N = params.ExtForces_N;
Ext_Moments_Nm = params.Extforces_N;
Mass_KG = params.Mass_KG;
BodyRates_radps = params.BodyRates_radps
LBE = params.LBE
Mass_kg = params.Mass_kg
g_mps = params.g_mps
GroundForcesAndMoments_N_Nm = params.GroundForcesAndMoments_N_Nm
AeroForcesAndMoments_N_Nm = params.AeroForcesAndMoments_N_Nm
ThurstForcesAndMoments_N_Nm = params.ThurstForcesAndMoments_N_Nm
InertiaTensor_kgm2 = params.InertiaTensor_kgm2;
V = params.V;
phi = EulerAngles_rad(1,1);
theta = EulerAngles_rad(2,1);
psi = EulerAngles_rad(3,1);


%Outuput
%Extforces_N
%Ext_Moments_Nm
%LoadFactors_N
%EulerRates_rad
%BodyVelocities_mps
%BodyRates_radps; EulerAngles_rad
%BodyRates_radps
%InertiaTensor_kgm2 = params.InertiaTensor_kgm2;

%Formulas

%InertiaTensor_kgm
 x = Ixx*p_dot - Ixz*r_dot - q*(Ixz*p - Izz*r) - Iyy*q*r;
 y = Iyy*q_dot + p*(Ixz*p - Izz*r) + r*(Ixx*p - Ixz*r);
 z = Izz*r_dot - Ixz*p_dot - q*(Ixx*p - Ixz*r) + Iyy*p*q;
iT = [x;y;z];

%AngularRotation
i = params.I;
j = params.J;
k = params.K;

Lpsi = plV.Lpsi;
Ltheta = plV.Ltheta;
Lphi = plV.Lphi;

pAngR = (Lpsi^-1)*(Ltheta^-1)*(Lphi^-1)*[i;j;k];

%BodyRates
pqr_dot = (InertiaTensor_kgm²^-1)*(InertiaTensor_kgm²- Ext_Moments_Nm);

pqr_radps = Integral(pqr_dot);

BodyRates_radps = (InertiaTensor_kgm²^-1)*(((InertiaTensor_kgm²*pqr_radps)x pqr_radps)- Ext_Moments_NmMoments_Nm));


u_dot = plV.Lpsi;
v_dot = plV.Ltheta;
w_dot = plV.Lphi;

pT = m*([u_dot; v_dot; w_dot] + cross([p;q;r],[u;v;w]));


%BodyVelocities
BodyVelocities = (ExtForces_N/Mass_KG) - (BodyRates_radps x uvw_mps)


%EulerRates
EulerRates_radps2 = BodyRates_radps + EulerAngles_rad;

%EulerRates_radps2 = EulerRates_radps2[((1;1)+(2;1)*(4;1)*(5;1))/(8;1))+(((3;1)*(7;1))*(5;1))/(8;1); ((2;1)*(7;1))-((3;1)*(4;1)); (((2;1)*(4;1))+((3;1)*(7;1)))/(8;1)];  (u(2)*u(4)+u(3)*u(7))/u(8)


%EulerAngles
EulerAngles_rad = integral(EulerRates_radps2);


%LBE(params)
Lpsi = [cos*EulerAngles_rad(psi) sin*EulerAngles_rad(psi) 0; -sin*EulerAngles_rad(psi) cos*EulerAngles_rad(psi) 0; 0 0 1];
Ltheta = [cos*EulerAngles_rad(theta) 0 -sin*EulerAngles_rad(theta); 0 1 0; sin*EulerAngles_rad(theta) 0 cos*EulerAngles_rad(2,1)];
Lphi = [1 0 0; 0 cos*EulerAngles_rad(phi) sin*EulerAngles_rad(phi); 0 -sin*EulerAngles_rad(phi) cos*EulerAngles_rad(phi)];

LBE = (Lphi*Ltheta*Lpsi)*.V;

%LBE = SICOS_Euler[(5,1)*(6,1) ((1,1)*(2,1)*(6,1))-((4,1)*(3,1)) ((6,1)*(2,1)*(4,1))+((1,1)*u(3,1)); (5,1)*(3,1) (1,1)*(2,1)*(3,1)+(4,1)*(6,1) ((4)*(2)*u(3))-((1)*(6)); -(2); (1,1)*(5,1) (4,1)*(5,1)];

W_N = [0;0;g_mps]*Mass_kg;

%ExtForces_N
Extforces_N = (LBE^T1_k * W_N) + (GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm);

%Ext_Moments_NmMoments
Ext_Moments_Nm = GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm;

%LoadFactors_N
LoadFactors_N = ((W_N^-1) * (GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm))*[1; 1; -1];

end
