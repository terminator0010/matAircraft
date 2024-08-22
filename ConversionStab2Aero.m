function [CxB, CyB, CzB, ClB, CmB, CnB, deltaXcg, Mcg] = ConversionStab2AeroCoeffi(params)
%Inputs
CL = params.CL;
CD = params.CD;
CY = params.CY;
CM = params.CM;
CN = params.CN;
M25 = params.M25;
X = params.X;
Y = params.Y;
CG = params.CG;

%Formula
CxB = CL*sin(Alpha) - CD*cos(Alpha);
CzB = -CL*cos(Alpha) - CD*sin(Alpha);
CyB = CY;

CmB = CM;
ClB = CL*cos(Alpha) - CN*sin(Alpha);
CnB = CL*sin(Alpha) + CN*cos(Alpha);

deltaXcg = CG-0,25;

%Mcg = M25 + deltaXcg*(-Z)-Zcg*(-X);
Mcg = M25 + X*Zcg - Z*deltaXcg;

end
