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




Alpha_radps = AlphaAngle_deg*(pi/180);
Beta_rapds = BetaAngle_deg*(pi/180);
Alphadot_radps;
TAS_mps;
rho_kgpm3;
BodyRates_radps;
Elevator_rapds = Elevator_deg*(pi/180);
Aileron_radps = Aileron_deg*(pi/180);
Rudder_radps = Rudder_deg*(pi/180);

CL_stabAxis = CL0+(CL_alpha*Alpha_radps)+(CL_elev*Elevator_rapds)+(CL_alphadot*Alphadot_radps)/(2*TAS_mps)+(q_radps*CL_q*c)/(2*TAS_mps);

CD_stabAxis = CD0+(|CD_beta|*|Alpha_radps|)+(|CD_elev|*|Elevator_rapds|);

CY_stabAxis = CY0+(CY_beta*beta_radps)+(CY_rud*Rudder_rapds)+(CY_Aileron*Aileron_radps)+(CY_r*r_radps*b)/(2*TAS_mps))+(p_radps*CY_p*b/(2*TAS_mps));

CI_stabAxis = CI0+(CI_beta*beta_radps)+(CI_rud*Rudder_rapds)+(CI_Aileron*Aileron_radps)+(CI_r*r_radps*b)/(2*TAS_mps))+(p_radps*CI_p*b/(2*TAS_mps));

CM_stabAxis = Cm0+(Cm_alpha*Alpha_radps)+(Cm_elev*Elevator_rapds)+(Cm_alphadot*Alphadot_radps/(2*TAS_mps))+(q_radps*Cm_q*c/(2*TAS_mps));

CN_stabAxis = CN0+(CN_beta*beta_radps)+(CN_rud*Rudder_rapds)+(CN_Aileron*Aileron_radps)+(CN_r*r_radps*b)/(2*TAS_mps))+(p_radps*CN_p*b/(2*TAS_mps));

C_stabAxis = [CL_stabAxis; CD_stabAxis; CY_stabAxis; CI_stabAxis; CM_stabAxis; CN_stabAxis];

end
