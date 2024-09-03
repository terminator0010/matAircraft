function C_stabAxis = Aerodynamics_Coefficient(params)

%Aerodynamic_Library_Coefficients
%Input
%TAS_mps
%AlphaAngle_deg
%BetaAngle_deg
%Alphadot_radps
%rho_kgpm3
%BodyRates_radps
%Elevator_deg
%Aileron_deg
%Rudder_deg


Alpha_radps = params.Alpha_angle_radps;
Beta_radps = params.Beta_angle_radps;
AlphaVelocities_radps = params.AlphaVelocities_radps;
TAS_mps = params.TAS_mps;
%perfGasEq = params.perfGasEq;
BodyRates_radps = params.BodyRates_radps;
Elevator_rapds = Elevator_deg*(pi/180);
Aileron_radps = Aileron_deg*(pi/180);
Rudder_radps = Rudder_deg*(pi/180);
p_radps = params.pqr_radps(1,1);
q_radps = params.pqr_radps(1,2);
r_radps = params.pqr_radps(1,3);


params.b = b;
params.c = c;
params.CL0 = CL0;
params.CL_alpha = CL_alpha;
params.CL_elev = CL_elev;
params.CL_AlphaDot = CL_AlphaDot;
params.CL_q = CL_q;
params.CD0 = CD0;
params.CD_alpha = CD_alpha;
params.CD_elev = CD_elev;
params.CY_beta = CY_beta;
params.CY_rud = CY_rud;
params.CY_ail = CY_Aileron;
params.CY_r = CY_r;
params.CY_p = CY_p;
params.CI_beta = CI_beta;
params.CI_rud = CI_rud;
params.CI_ail = CI_Aileron;
params.CI_r = CI_r;
params.CI_p = CI_p;
params.Cm0 = Cm0;
params.Cm_alpha = Cm_alpha;
params.Cm_elev = Cm_elev;
params.Cm_AlphaDot = Cm_AlphaDot;
params.Cm_q = Cm_q;
params.Cn_beta = CN_beta;
params.Cn_rud = CN_rud;
params.Cn_ail = CN_Aileron;
params.Cn_r = CN_r;
params.Cn_p = CN_p;



CL_stabAxis = CL0+(CL_alpha*Alpha_radps)+(CL_elev*Elevator_rapds)+(CL_AlphaDot*AlphaVelocities_radps)/(2*TAS_mps)+(q_radps*CL_q*c)/(2*TAS_mps);

CD_stabAxis = CD0+(abs(CD_alpha)*abs(Alpha_radps))+(abs(CD_elev)*abs(Elevator_rapds));

CY_stabAxis = CY0+(CY_beta*Beta_radps)+(CY_rud*Rudder_radps)+(CY_Aileron*Aileron_radps)+(CY_r*r_radps*b)/(2*TAS_mps)+(p_radps*CY_p*b)/(2*TAS_mps);

CI_stabAxis = CI0+(CI_beta*Beta_radps)+(CI_rud*Rudder_radps)+(CI_Aileron*Aileron_radps)+(CI_r*r_radps*b)/(2*TAS_mps)+(p_radps*CI_p*b)/(2*TAS_mps);

CM_stabAxis = Cm0+(Cm_alpha*Alpha_radps)+(Cm_elev*Elevator_rapds)+(Cm_AlphaDot*AlphaVelocities_radps)/(2*TAS_mps)+(q_radps*Cm_q*c)/(2*TAS_mps);

CN_stabAxis = CN0+(CN_beta*Beta_radps)+(CN_rud*Rudder_radps)+(CN_Aileron*Aileron_radps)+(CN_r*r_radps*b)/(2*TAS_mps)+(p_radps*CN_p*b)/(2*TAS_mps);

C_stabAxis = [CL_stabAxis; CD_stabAxis; CY_stabAxis; CI_stabAxis; CM_stabAxis; CN_stabAxis];

end
