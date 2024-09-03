function AeroForcesMoments = AeroforcesandMoments(params)

%Input
perfGasEq = params.perfGasEq;
TAS_mps = params.TAS_mps;
WingArea_m2 = params.S;
CG_mac = params.CG_Mac;
Zcg_m = params.Zcg_m;
Ycg_m = params.Ycg_m;
N25 = params.N25;
L25 = params.L25;
b = params.b;
c = params.c;
deltaXcg = params.deltaXcg;
CxB = params.CxB;
CyB = params.CyB;
CzB = params.CzB;
ClB = params.ClB;
CmB = params.CmB;
CnB = params.CnB;
Mcg = params.Mcg;
Ixx = params.Ixx;
Iyy = params.Iyy;
Izz = params.Izz;


%Yaw Cg
Ncg = N25 + Iyy*deltaXcg + Ixx*Ycg_m

%Rolling Cg
Lcg = L25 + Iyy*Zcg_m - Izz*Ycg_m

dynP_Pa = (perfGasEq/2) * TAS_mps^2;

Qstab = dynP_Pa * WingArea_m2;

X_Aero_N = Qstab * CxB;
Y_Aero_N = Qstab * CyB;
Z_Aero_N = Qstab * CzB;

Fcn = (CG_mac-0.25)*c;

M_Aero_N = ((CmB*Qstab)*c) - (Z_Aero_N * Fcn) + (X_Aero_N * Zcg_m);

N_Aero_N = (Y_Aero_N * Fcn) + (X_Aero_N * Ycg_m) + ((CnB*Qstab)*b);

L_Aero_N = ((ClB*Qstab)*b) - (Z_Aero_N * Ycg_m) - (Y_Aero_N * Zcg_m);

AeroForcesMoments = [X_Aero_N; Y_Aero_N; Z_Aero_N; M_Aero_N; N_Aero_N; L_Aero_N];

end
