function [Extforces_N, Ext_Moments_Nm, LoadFactors_N] = extForcesAndMoments (params)
%Input
LBE = params.LBE
Mass_kg = params.Mass_kg
g_mps = params.g_mps
GroundForcesAndMoments_N_Nm = params.GroundForcesAndMoments_N_Nm
AeroForcesAndMoments_N_Nm = params.AeroForcesAndMoments_N_Nm
ThurstForcesAndMoments_N_Nm = params.ThurstForcesAndMoments_N_Nm

%Outuput
%Extforces_N
%Ext_Moments_Nm
%LoadFactors_N


%Formulas
W_N = [0;0;g_mps]*Mass_kg;
Extforces_N = (LBE^T1_k * W_N) + (GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm);
Ext_Moments_Nm = GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm;
LoadFactors_N = ((W_N^-1) * (GroundForcesAndMoments_N_Nm + AeroForcesAndMoments_N_Nm  + ThurstForcesAndMoments_N_Nm))*[1; 1; -1];

end



