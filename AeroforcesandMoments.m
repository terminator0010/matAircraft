function [X_Aero_N, Y_Aero_N, Z_Aero_N, M_Aero_N, N_Aero_N, L_Aero_N, dynP_Pa] = AeroforcesandMoments(params)

%Input
perfGasEq = params.perfGasEq;
TAS_mps = params.TAS_mps;
WingArea_m2 = params.WingArea_m2;
CG_Mac = params.CG_Mac;
Zcg_m = params.Zcg_m;
Ycg_m = params.Ycg_m;
b = params.b;
c = params.c;
CxB = params.CxB;
CyB = params.CyB;
CzB = params.CzB;
ClB = params.ClB;
CmB = params.CmB;
CnB = params.CnB;
Mcg = params.Mcg;
X = params.PositionInertial_m(1,1);
Y = params.PositionInertial_m(2,1);
Z = params.PositionInertial_m(3,1);
deltaXcg = params.deltaXcg;




%Yaw Cg

%M25 = ((CzB*dynP_Pa)*deltaXcg)+(Zcg_m*CxB)+CG_Mac;


%Rolling Cg
%Lcg = L25 + Y*Zcg_m - Z*Ycg_m

dynP_Pa = (perfGasEq/2) * TAS_mps^2;

Qstab = dynP_Pa * WingArea_m2;

X_Aero_N = Qstab * CxB;
Y_Aero_N = CyB;
Z_Aero_N = Qstab * CzB;

Fcn = (CG_Mac-0.25)*c;

M_Aero_N = ((CmB*Qstab)*c) - (Z_Aero_N * Fcn) + (X_Aero_N * Zcg_m);

N_Aero_N = (Y_Aero_N * Fcn) + (X_Aero_N * Ycg_m) + ((CnB*Qstab)*b);

L_Aero_N = ((ClB*Qstab)*b) - (Z_Aero_N * Ycg_m) - (Y_Aero_N * Zcg_m);

end
