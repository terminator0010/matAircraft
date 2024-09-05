%% INPUT DATA
function [mass, cg, Alt_ft, KCAS, Gamma_deg, Phi_deg, Throttle] = trim_input(params)
mass = 18000;
cg = 0.47;
Alt_ft = 40000;
KCAS = 230;
Gamma_deg = 0;
Phi_deg = 0;
Throttle = -1; % if throtle is -1, then will float throttle
%{
%% Search for a specified operating point for the model - ACFT.
%
% This MATLAB script is the command line equivalent of the trim model
% tab in linear analysis tool with current specifications and options.
% It produces the exact same operating points as hitting the Trim button.

% MATLAB(R) file generated by MATLAB(R) 9.0 and Simulink Control Design (TM) 4.3.
%
% Generated on: 25-Jul-2020 14:06:41

%% Specify the model name
model = 'ACFT';

%% Create the operating point specification object.
opspec = operspec(model);

%% Set the constraints on the states in the model.
% - The defaults for all states are Known = false, SteadyState = true,
%   Min = -Inf, and Max = Inf.

% State (1) - ACFT/EqM/BodyRates/p_radps
% - Default model initial conditions are used to initialize optimization.

% State (2) - ACFT/EqM/BodyRates/q_radps
% - Default model initial conditions are used to initialize optimization.

% State (3) - ACFT/EqM/BodyRates/r_radps
% - Default model initial conditions are used to initialize optimization.

% State (4) - ACFT/EqM/BodyVelocities/u_mps
opspec.States(4).x = KCAS/2;

% State (5) - ACFT/EqM/BodyVelocities/v_mps
% - Default model initial conditions are used to initialize optimization.

% State (6) - ACFT/EqM/BodyVelocities/w_mps
% - Default model initial conditions are used to initialize optimization.

% State (7) - ACFT/EqM/EulerAngles/PHI_rad
% - Default model initial conditions are used to initialize optimization.

% State (8) - ACFT/EqM/EulerAngles/PSI_rad
% - Default model initial conditions are used to initialize optimization.
if Phi_deg ~= 0
    opspec.States(8).SteadyState = false;
else
    opspec.States(8).SteadyState = true;
end

% State (9) - ACFT/EqM/EulerAngles/THETA_rad
% - Default model initial conditions are used to initialize optimization.

% State (10) - ACFT/EqM/InertialData/XI_m
% - Default model initial conditions are used to initialize optimization.
opspec.States(10).SteadyState = false;

% State (11) - ACFT/EqM/InertialData/YI_m1
% - Default model initial conditions are used to initialize optimization.
opspec.States(11).SteadyState = false;

% State (12) - ACFT/EqM/InertialData/ZI_m
% - Default model initial conditions are used to initialize optimization.
opspec.States(12).SteadyState = false;

%% Set the constraints on the inputs in the model.
% - The defaults for all inputs are Known = false, Min = -Inf, and
% Max = Inf.

% Input (1) - ACFT/Mass_kg
opspec.Inputs(1).u = mass;
opspec.Inputs(1).Known = true;

% Input (2) - ACFT/CG_mac
opspec.Inputs(2).u = cg;
opspec.Inputs(2).Known = true;

% Input (3) - ACFT/Elevator_deg
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(3).Min = -30;
opspec.Inputs(3).Max = 30;

% Input (4) - ACFT/Aileron_deg
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(4).Min = -40;
opspec.Inputs(4).Max = 40;

% Input (5) - ACFT/Rudder_deg
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(5).Min = -30;
opspec.Inputs(5).Max = 30;

% Input (6) - ACFT/Throttle
if Throttle~=-1;
opspec.Inputs(6).u = Throttle;
opspec.Inputs(6).Min = 0;
opspec.Inputs(6).Max = 1;
opspec.Inputs(6).Known = true;
else
opspec.Inputs(6).Min = 0;
opspec.Inputs(6).Max = 1;
opspec.Inputs(6).Known = false;
end

% Input (7) - ACFT/WindX
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(7).Known = true;

% Input (8) - ACFT/WindY
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(8).Known = true;

% Input (9) - ACFT/WindZ
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(9).Known = true;

%% Set the constraints on the outputs in the model.
% - The defaults for all outputs are Known = false, Min = -Inf, and
% Max = Inf.

% Output (1) - ACFT/KCAS
opspec.Outputs(1).y = KCAS;
opspec.Outputs(1).Known = true;

% Output (2) - ACFT/Mach
% - Default model initial conditions are used to initialize optimization.

% Output (3) - ACFT/KTAS
% - Default model initial conditions are used to initialize optimization.

% Output (4) - ACFT/Beta_deg
% - Default model initial conditions are used to initialize optimization.
opspec.Outputs(4).Known = true;

% Output (5) - ACFT/Alpha_deg
% - Default model initial conditions are used to initialize optimization.

% Output (6) - ACFT/AlphaDot_radps
% - Default model initial conditions are used to initialize optimization.

% Output (7) - ACFT/u_mps
% - Default model initial conditions are used to initialize optimization.

% Output (8) - ACFT/v_mps
% - Default model initial conditions are used to initialize optimization.

% Output (9) - ACFT/w_mps
% - Default model initial conditions are used to initialize optimization.

% Output (10) - ACFT/p_radps
% - Default model initial conditions are used to initialize optimization.

% Output (11) - ACFT/q_radps
% - Default model initial conditions are used to initialize optimization.

% Output (12) - ACFT/r_radps
% - Default model initial conditions are used to initialize optimization.

% Output (13) - ACFT/Phi_deg
opspec.Outputs(13).y = Phi_deg;
opspec.Outputs(13).Known = true;

% Output (14) - ACFT/Theta_deg
% - Default model initial conditions are used to initialize optimization.

% Output (15) - ACFT/Psi_deg
% - Default model initial conditions are used to initialize optimization.

% Output (16) - ACFT/X_m
% - Default model initial conditions are used to initialize optimization.

% Output (17) - ACFT/Y_m
% - Default model initial conditions are used to initialize optimization.

% Output (18) - ACFT/PresAlt_ft
opspec.Outputs(18).y = Alt_ft;
opspec.Outputs(18).Known = true;

% Output (19) - ACFT/GSpeed_kt
% - Default model initial conditions are used to initialize optimization.

% Output (20) - ACFT/Gamma_deg
% - Default model initial conditions are used to initialize optimization.
if Throttle~=-1;
opspec.Outputs(20).Known = false;
else
    opspec.Outputs(20).Known = true;
    opspec.Outputs(20).y = Gamma_deg;
end


% Output (21) - ACFT/Track_deg
% - Default model initial conditions are used to initialize optimization.

% Output (22) - ACFT/nx
% - Default model initial conditions are used to initialize optimization.

% Output (23) - ACFT/ny
% - Default model initial conditions are used to initialize optimization.

% Output (24) - ACFT/nz
% - Default model initial conditions are used to initialize optimization.

% Output (25) - ACFT/Thrust_N
% - Default model initial conditions are used to initialize optimization.


%% Create the options
opt = findopOptions('DisplayReport','iter');

%% Perform the operating point search.
[op,opreport] = findop(model,opspec,opt);
%}

end
