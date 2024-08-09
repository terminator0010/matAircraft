%function [EulerAngles_rad, BodyVelocities, Extforces_N, Ext_Moments_Nm, LoadFactors_N, LBE, PositionInertial_m, GroundSpeed_mps, Gamma_angle_deg, Track_angle_deg] = eQMotion(params)

function [InertiaTensor_kgm2, pqr_dot, pqr_radps, BodyRates_radps, I, BodyVelocities, EulerRates_radps2, EulerAngles_rad, LBE, uvw_dot, uvw_mps, pT] = eQMotion(params)

%function [InertiaTensor_kgm2, pqr_dot, pqr_radps, I]  = eQMotion(params)

%Input
ExtForces_N = params.ExtForces_N;
Ext_Moments_Nm = params.Ext_Moments_Nm;
Mass_KG = params.Mass_KG;
BodyRates_radps = params.BodyRates_radps;
g_mps = params.g_mps;
GroundForcesAndMoments_N_Nm = params.GroundForcesAndMoments_N_Nm;
AeroForcesAndMoments_N_Nm = params.AeroForcesAndMoments_N_Nm;
ThurstForcesAndMoments_N_Nm = params.ThurstForcesAndMoments_N_Nm;
InertiaTensor_kgm2 = params.InertiaTensor_kgm2;
V = params.V;
Rotation = params.Rotation;
%EulerAngles_rad = [phi = params.EulerAngles_rad(1,1); theta = params.EulerAngles_rad(2,1); psi = params.EulerAngles_rad(3,1)];
phi = params.EulerAngles_rad(1,1);
theta = params.EulerAngles_rad(2,1);
psi = params.EulerAngles_rad(3,1);
EulerAngles_rad = [phi; theta; psi];
InertiaTensor_kgm2 = params.InertiaTensor_kgm2;
pqr_radps = params.pqr_radps;
pqr_dot = params.pqr_dot;
uvw_mps = params.uvw_mps;
t1_k = params.t1_k;


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


 %x = Ixx*p_dot - Ixz*r_dot - q*(Ixz*p - Izz*r) - Iyy*q*r;
 %y = Iyy*q_dot + p*(Ixz*p - Izz*r) + r*(Ixx*p - Ixz*r);
 %z = Izz*r_dot - Ixz*p_dot - q*(Ixx*p - Ixz*r) + Iyy*p*q;

 InertiaTensor_kgm2 = [InertiaTensor_kgm2(1,1)*pqr_dot(1,1) - InertiaTensor_kgm2(1,3)*pqr_dot(2,1) - pqr_radps(3,1)*(InertiaTensor_kgm2(1,3)*pqr_radps(1,1) - InertiaTensor_kgm2(3,3)*pqr_radps(2,1)) - InertiaTensor_kgm2(2,2)*pqr_radps(3,1)*pqr_radps(1,1);
                       InertiaTensor_kgm2(2,2)*pqr_dot(2,1) + pqr_radps(1,1)*(InertiaTensor_kgm2(1,3)*pqr_radps(1,1) - InertiaTensor_kgm2(3,3)*pqr_radps(3,1)) + pqr_radps(3,1)*(InertiaTensor_kgm2(1,1)*pqr_radps(1,1)) - (InertiaTensor_kgm2(1,3)*pqr_radps(3,1));
                       InertiaTensor_kgm2(3,3)*pqr_dot(3,1) - InertiaTensor_kgm2(1,3)*pqr_dot(1,1) - pqr_radps(2,1)*(InertiaTensor_kgm2(1,1)* pqr_radps(1,1) - InertiaTensor_kgm2(1,3)*pqr_radps(3,1)) + InertiaTensor_kgm2(2,2)*pqr_radps(1,1)*pqr_radps(2,1)];



%BodyRates
pqr_dot = (InertiaTensor_kgm2.*-1).*(InertiaTensor_kgm2 - Ext_Moments_Nm);

pqr_radps = trapz(pqr_dot,2);

I = InertiaTensor_kgm2.*pqr_radps;

%I = cross(pqr_radps,InertiaTensor_kgm2*pqr_radps);

%I = InertiaTensor_kgm2.*-1;

BodyRates_radps = (InertiaTensor_kgm2.*-1).*(cross(pqr_radps,I))- Ext_Moments_Nm;





%BodyVelocities
BodyVelocities = (ExtForces_N/Mass_KG) - (cross(BodyRates_radps, uvw_mps));


%EulerRates
EulerRates_radps2 = BodyRates_radps + EulerAngles_rad;

%EulerRates_radps2 = EulerRates_radps2[((1;1)+(2;1)*(4;1)*(5;1))/(8;1))+(((3;1)*(7;1))*(5;1))/(8;1); ((2;1)*(7;1))-((3;1)*(4;1)); (((2;1)*(4;1))+((3;1)*(7;1)))/(8;1)];  (u(2)*u(4)+u(3)*u(7))/u(8);


%EulerAngles


EulerAngles_rad = trapz(EulerRates_radps2, 2);


%LBE(params)
Lpsi = [cos(EulerAngles_rad(psi)) sin(EulerAngles_rad(psi)) 0; -sin(EulerAngles_rad(psi)) cos(EulerAngles_rad(psi)) 0; 0 0 1];
Ltheta = [cos(EulerAngles_rad(theta)) 0 -sin(EulerAngles_rad(theta)); 0 1 0; sin(EulerAngles_rad(theta)) 0 cos(EulerAngles_rad(2,1))];
Lphi = [1 0 0; 0 cos(EulerAngles_rad(phi)) sin(EulerAngles_rad(phi)); 0 -sin(EulerAngles_rad(phi)) cos(EulerAngles_rad(phi))];

LBE = (Lphi*Ltheta*Lpsi).*V;

%LBE = SICOS_Euler[(5,1)*(6,1) ((1,1)*(2,1)*(6,1))-((4,1)*(3,1)) ((6,1)*(2,1)*(4,1))+((1,1)*u(3,1)); (5,1)*(3,1) (1,1)*(2,1)*(3,1)+(4,1)*(6,1) ((4)*(2)*u(3))-((1)*(6)); -(2); (1,1)*(5,1) (4,1)*(5,1)];


uvw_dot = [Lphi(1,1); Ltheta(2,2); Lpsi(3,3)];
uvw_mps = trapz(uvw_dot, 2);

pT = Mass_KG*(uvw_dot + cross(pqr_radps, uvw_mps));


%Angular Rotation
pAngR = (Lpsi^-1)*(Ltheta^-1)*(Lphi^-1).*Rotation;

W_N = [0;0;g_mps]*Mass_KG;
%{
%ExtForces_N
Extforces_N = (LBE^T1_k * W_N) + (GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm);


%Ext_Moments_NmMoments
Ext_Moments_Nm = GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm;

%LoadFactors_N
LoadFactors_N = ((W_N^-1) * (GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm))*[1; 1; -1];

%Inertial Data
InertialVelocity_mps = LBE * BodyVelocities;
PositionInertial_m = trapz(InertialVelocity_mps,2);


%Trajectory Data
%Input
%Ve_dot (Xe_dot; Ye_dot; Ze_dot)
%InertialVelocity_mps


%Output
%Gamma_angle_deg
%Track_angle_deg
%GroundSpeed_mps

%Formula
Xe_dot = InertialVelocity_mps(1,1);
Ye_dot = InertialVelocity_mps(1,2);
Ze_dot = InertialVelocity_mps(1,3);

%GroundSpeed_mps
GroundSpeed_mps = sqrt(Xe_dot^2+Ye_dot^2);

%Gamma_Angle
%gamma_angle = (tg^(-1))*((-1*Ze_dot)/(sqrt(Xe_dot^2+Ye_dot^2)))
Gamma_angle_deg = (atan2((-1*Ze_dot)/GroundSpeed_mps))*(180*pi);

%Track_angle
Track_angle_deg = atan2(Ye_dot/Xe_dot);
%}

end

