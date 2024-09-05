function [CxB, CyB, CzB, ClB, CmB, CnB, deltaXcg, Mcg] = ConversionStab2Aero(params)
%Inputs
CL0 = params.CL0;
CD0 = params.CD0;
CY0 = params.CY0;
CM = params.CM;
CN = params.CN;
M25 = params.M25;
X = params.PositionInertial_m(1,1);
Y = params.PositionInertial_m(2,1);
Z = params.PositionInertial_m(3,1);
CG = params.CG;

%Formula
CxB = CL0*sin(Alpha) - CD0*cos(Alpha);
CzB = -CL0*cos(Alpha) - CD0*sin(Alpha);
CyB = CY0;

CmB = CM;
ClB = CL*cos(Alpha) - CN*sin(Alpha);
CnB = CL*sin(Alpha) + CN*cos(Alpha);

deltaXcg = CG-0,25;

%Mcg = M25 + deltaXcg*(-Z)-Zcg*(-X);
Mcg = M25 + X*Zcg - Z*deltaXcg;

end
