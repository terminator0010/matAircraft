function BodyRates

Input
Ext_Moments_Nm
InertiaTensor_kgm²

Output
BodyRates_radps

Formulas
BodyRates_radps = ([InertialTensor_kgm²]^-1)*((([InertialTensor_kgm²]*[p_radps; q_radps; r_radps])x[p_radps; q_radps; r_radps])- Moments_Nm))
[p_dot; q_dot; r_dot] = ([InertialTensor_kgm²]^-1)*([InertialTensor_kgm²]- Moments_Nm)
[p_radps; q_radps; r_radps] = Integral[p_dot; q_dot; r_dot]
end

