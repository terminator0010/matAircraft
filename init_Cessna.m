% Cessna 172
% Mass properties
Ixx = 1285.3;                                                               % kg.m2
Iyy = 1824.9;                                                               % kg.m2
Izz = 2666.9;                                                               % kg.m2
Ixz = 0;                                                                    % kg.m2
Inertia = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];                                % Inertia tensor
y_cg = 0;                                                                   % m (+ is to the right of ref. plane)
z_cg = 0.2;                                                                 % m (+ is below the ref. plane)
g = 9.80665;                                                                % m/s2

% Engine Data
V_reference_mps = 51.4;                                                     % m/s
rho_reference_kgpm3 = 1.225;                                                % kg/m3
nv = -1;
nrho = 0.75;
alfa_deg = 1.0;                                                            % Thrust angle - deg
xf_m = 1.0;                                                                 % Enigne distance on x axis from approximate CG position - m (+ is FWD of the CG)
zf_m =  0.0;                                                                % Enigne distance on z axis from approximate CG position - m (+ is below the ref. plane)
Tmax = 2070;                                                                % Thrust - N

% Aircraft Data
S = 16.2;                                                                   % Wing Area - m2
c = 4.9*0.3048;                                                             % Mean Aerodynamic Chord - m
b = 36*0.3048;                                                              % Wing Span - m

% Aerodynamic Data

% Lift Coefficient Derivatives
CL0 = 0.307;
CL_alpha = 4.41;
CL_elev = 0.43;
CL_AlphaDot = 0.0;
CL_q = 3.9;
% Drag Coefficient Derivatives
CD0 = 0.0270;
CD_alpha = 0.26;
CD_elev = 0.0;
% Side force Coefficient Derivatives
CY_beta = -0.404;
CY_rud = 0.187;
CY_ail = 0;
CY_r = 0.214;
CY_p = -0.075;
% Rolling moment Derivatives
Cl_beta = -0.0923;
Cl_rud = 0.0147;
Cl_ail = -0.229;
Cl_r = 0.0798;
Cl_p = -0.484;
% Pitching moment Derivatives
Cm0 = 0.04;
Cm_alpha = -0.613;
Cm_elev = -1.122;
Cm_AlphaDot = -7.27;
Cm_q = -12.4;
% Yawing moment Derivatives
Cn_beta = 0.0587;
Cn_rud = -0.0645;
Cn_ail = -0.0216;
Cn_r = -0.0937;
Cn_p = -0.0278;

%Trimming
Trimming=0;
