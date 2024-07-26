function BodyRates_radps = BodyRates(params)

%Input
Ext_Moments_Nm = params.Extforces_N;
InertiaTensor_kgm2 = params.InertiaTensor_kgm2;

%Output
BodyRates_radps

%Formulas
BodyRates_radps = ([InertiaTensor_kgm²]^-1)*((([InertiaTensor_kgm²]*[p_radps; q_radps; r_radps])x[p_radps; q_radps; r_radps])- Moments_Nm));

[p_dot; q_dot; r_dot] = ([InertiaTensor_kgm²]^-1)*([InertiaTensor_kgm²]- Moments_Nm);

[p_radps; q_radps; r_radps] = Integral[p_dot; q_dot; r_dot];

end

