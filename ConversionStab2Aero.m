function [CxB, CyB, CzB, ClB, CmB, CnB, deltaXcg] = ConversionStab2Aero(params)
%Inputs
CL_stabAxis = params.CL_stabAxis;
CD_stabAxis = params.CD_stabAxis;
CY_stabAxis = params.CY_stabAxis;
CN_stabAxis = params.CN_stabAxis;
CM_stabAxis = params.CN_stabAxis;
Alpha_angle_radps = params.Alpha_angle_radps;
X = params.PositionInertial_m(1,1);
Y = params.PositionInertial_m(2,1);
Z = params.PositionInertial_m(3,1);
CG_Mac = params.CG_Mac;
Zcg_m = params.Zcg_m;
Ycg_m = params.Ycg_m;
c = params.c;
b = params.b;
dynP_Pa = params.dynP_Pa;

%Formula
CxB = CL_stabAxis*sin(Alpha_angle_radps) - CD_stabAxis*cos(Alpha_angle_radps);
CzB = -CL_stabAxis*cos(Alpha_angle_radps) - CD_stabAxis*sin(Alpha_angle_radps);
CyB = CY_stabAxis;

CmB = CM_stabAxis;
ClB = CL_stabAxis*cos(Alpha_angle_radps) - CN_stabAxis*sin(Alpha_angle_radps);
CnB = CL_stabAxis*sin(Alpha_angle_radps) + CN_stabAxis*cos(Alpha_angle_radps);

c = dynP_Pa*CmB;
b = dynP_Pa*CnB;

deltaXcg = ((CG_Mac-0.25)*c);

%Mcg = M25 + deltaXcg*(-Z)-Zcg*(-X);
%M25 = ((CzB*dynP_Pa)*deltaXcg)+(Zcg_m*CxB)+CG_Mac;
%N25 = ((CyB*dynP_Pa)*deltaXcg)+(Ycg_m*CxB)+CG_Mac;

end
