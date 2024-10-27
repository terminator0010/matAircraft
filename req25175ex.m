
%% SIMULATION
% Writing initial inputs vector.
inputs = getinputstruct(opreport);                                          % Reads the model inputs
utin = zeros(size(inputs.signals,2),1);                                     % initialize with zeros the inputs vector - this is optional but reduces computational time of the code.
for i = 1:size(inputs.signals,2)
    utin(i,:) = inputs.signals(i).values;                                   % Creates a vector with the trimmed values for the inputs.
end

% Total simulation time - [s];
TF=1500;
% Time data for the input;
t = [0	5	55	105	155	205	255	305	355	405	455	505	555	605	655	705	850	1500];    % This vector can be any kind of time vector to create the input.

% Create an input vector of the same size as t;
%First line is the time stamp itself (simulink coding!)
ut = zeros(size(t,2),size(inputs.signals,2));                               % initialize with zeros the inputs vector - this is optional but reducis computational time of the code.
for i=1:size(t,2)
    ut(i,1) = t(i);                                                         % Simulink default. The first line is the time vector.
    for j=1:size(inputs.signals,2)
        ut(i,j+1) = utin(j);                                                % Initial input vectors with the same size as t.
    end
end

% ELEVATOR DOUBLET
                                                                            %amplitude [deg]
%Control Surface index
for i=1:size(opreport.Inputs,2)
    intr{i,1} = opreport.Inputs(i).Block;                                   % Create a structure with the names of the inputs
end
in_ind = [1:size(opreport.Inputs,2)];                                       % Create a vector with the indices of the inputs.
ind_control = in_ind(logical(strcmp('ACFT/Elevator_deg',intr)));            % Finds the position of the Elevator input.

% input command with the same size as t; You can draw anything here.
ut(:,ind_control+1) = [ut(1,ind_control+1) ut(1,ind_control+1) (ut(1,ind_control+1)-0.24) (ut(1,ind_control+1)-0.24) (ut(1,ind_control+1)-0.594) (ut(1,ind_control+1)-0.594) (ut(1,ind_control+1)-1.148) (ut(1,ind_control+1)-1.148) (ut(1,ind_control+1)-1.55) (ut(1,ind_control+1)-1.55) (ut(1,ind_control+1)-2.1) (ut(1,ind_control+1)-2.1) (ut(1,ind_control+1)-2.85) (ut(1,ind_control+1)-2.85) (ut(1,ind_control+1)-3.94) (ut(1,ind_control+1)-3.94) (ut(1,ind_control+1)) (ut(1,ind_control+1))]' ;

%Simulation command
% [tout,xout,yout]=sim('ACFT',TF,simset('InitialState',getstatestruct(opreport),'Solver','ode4','FixedStep',1),ut);

simOutput=sim('ACFT',TF,simset('InitialState',getstatestruct(opreport),'Solver','ode4','FixedStep',0.01),ut);

% Access specific results from the 'simOutput' variable
tout = simOutput.tout; % Time vector
yout_dataset = simOutput.yout; % Output variables as a 1x1 dataset


% Extract individual signals and store as a matrix (assuming 25 output signals)
  num_signals = 25;
  yout = zeros(length(tout), num_signals); % Initialize the matrix
  
  for i = 1:num_signals
      signal_data = yout_dataset.getElement(i).Values.Data; % Get the data of the i-th signal
      yout(:, i) = signal_data;
  end

%Plot
%Read Outputs
for i=1:size(opreport.Outputs,1)
    outr{i,1} = opreport.Outputs(i).Block;                                   % Create a structure with the names of the outputs.
end
out_ind = [1:size(opreport.Outputs,1)];                                      % Create a vector with the indices of the outputs.

grid;
figure(1);
subplot(421);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/Alpha_deg',outr))))); xlabel('Time - [s]');ylabel('Alpha - [deg]')
subplot(422);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/KCAS',outr))))); xlabel('Time - [s]');ylabel('KCAS')
subplot(423);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/Theta_deg',outr))))); xlabel('Time - [s]');ylabel('Theta - [deg]')
subplot(424);plot(t,ut(:,ind_control+1)); xlabel('Time - [s]');ylabel('Elevator - [deg]')
subplot(425);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/q_degps',outr))))); xlabel('Time - [s]');ylabel('Pitch Rate - [deg/s]')
subplot(426);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/Thrust_N',outr))))); xlabel('Time - [s]');ylabel('Thrust - [N]')
subplot(427);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/PresAlt_ft',outr))))); xlabel('Time - [s]');ylabel('Altitude - [ft]')
subplot(428);plot(tout,yout(:,out_ind(logical(strcmp('ACFT/nz',outr))))); xlabel('Time - [s]');ylabel('Load Factor (Nz) - [g]')
