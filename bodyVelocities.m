function BodyVelocities

Input
ExtForces_N
Mass_KG (m)
BodyRates_radps

Output
BodyVelocities_mps

Formulas
BodyVelocities = (ExtForces_N/Mass_KG) - (BodyRates_radps x [u_mps; v_mps; v_mps])

[u_dot; v_dot; q_dot] = (ExtForces_N/Mass_KG) - (BodyRates_radps)

[u_mps; v_mps; v_mps] = Integral[u_dot; v_dot; q_dot]


EulerAngles
Input
BodyRates_radps
EulerAngles_rad

Output
EulerRates_rad
LBE

Formulas
EulerRates_radps2 = [BodyRates_radps] + [EulerAngles_rad]

EulerRates_radps2 = EulerRates_radps2[(1;1)+((2;1)*(4;1))*((5;1)/(8;1))+((3;1)*(7;1))*((5;1)/(8;1)); ((2;1)*(7;1))-((3;1)*(4;1)); (((2;1)*(4;1))+((3;1)*(7;1))/(8;1)];  (u(2)*u(4)+u(3)*u(7))/u(8)

EulerAngles_rad = integral(EulerRates_radps2);

LBE = EulerAngles_rad[u(5)*u(6) (u(1)*u(2)*u(6))-(u(4)*u(3)) (u(6)*u(2)*u(4))+(u(1)*u(3)); u(5)*u(3) u(1)*u(2)*u(3)+u(4)*u(6) (u(4)*u(2)*u(3))-(u(1)*u(6)); -u(2); u(1)*u(5) u(4)*u(5)]

end
