function [Extforces_N; Ext_Moments_Nm; LoadFactors_N] = extForcesAndMoments (params)
LEB = [x; y; z];

m = params.m;
g = params.g;

m_N = [0;0;m];
Mass_KG = [m_N*.g;]


GroundForcesAndMoments_N_Nm
AeroForcesAndMoments_N_Nm
ThrustForcesAndMoments_N_Nm


Extforces_N = ([LBE^T1_k] * [Mass_kg]) + ([GroundForcesAndMoments_N_Nm] + [AeroForcesAndMoments_N_Nm]  + [ThurstForcesAndMoments_N_Nm])
Ext_Moments_Nm = ([GroundForcesAndMoments_N_Nm] + [AeroForcesAndMoments_N_Nm]  + [ThurstForcesAndMoments_N_Nm])
LoadFactors_N = ((([Mass_kg]^-1) * ([GroundForcesAndMoments_N_Nm] + [AeroForcesAndMoments_N_Nm]  + [ThurstForcesAndMoments_N_Nm]))* [1; 1; -1])

end



