function [CxB, CyB, CzB, ClB, CmB, CnB, deltaXcg, Mcg] = ConversionStab2Aero(params)
%Inputs
CL0 = params.CL0;
CD0 = params.CD0;
C_stabAxis = params.C_stabAxis;
CM = params.CM;
CN = params.CN;
M25 = params.M25;
params.Alpha_angle_radps
X = params.PositionInertial_m(1,1);
Y = params.PositionInertial_m(2,1);
Z = params.PositionInertial_m(3,1);
CG = params.CG;

%Formula
CxB = CL0*sin(Alpha_angle_radps) - CD0*cos(Alpha_angle_radps);
CzB = -CL0*cos(Alpha_angle_radps) - CD0*sin(Alpha_angle_radps);
CyB = C_stabAxis(3,1);

CmB = CM;
ClB = CL*cos(Alpha_angle_radps) - CN*sin(Alpha_angle_radps);
CnB = CL*sin(Alpha_angle_radps) + CN*cos(Alpha_angle_radps);

deltaXcg = CG-0,25;

%Mcg = M25 + deltaXcg*(-Z)-Zcg*(-X);
Mcg = M25 + X*Zcg - Z*deltaXcg;

end
