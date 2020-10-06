initModel5MWNRELBianchi;
%addpath(fullfile(fileparts(pwd),'datain'));
parentDir = fileparts(fileparts(fileparts(pwd)));
addpath(genpath(parentDir))
%addpath(fullfile(fileparts(parentDir),'datain'));
load('workSpace.mat')
load('OutDataSweep.mat');
OutTable = array2table(OutData,'VariableNames', OutList);
idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

%%
% ny = 2; nx = 4; nu=2;
% Carray1 = [eye(ny),zeros(ny,nx)];
% Darray1 = zeros(ny,nu);
% figure;
% A1 =  A(:,:,1);
% B1 =  B(:,:,1);
% SYSD1 = ss(A1,B1,Carray1,Darray1,DT);
% impulse(SYSD1,'k--');

%%
% Plot input
time = OutTable.Time(1: length(uDot));
figure; 
subplot(4,1,1); plot(time, uDot(:,1)); axis tight; grid on; 
ylabel('betaDot [rad/s]')
title('Simulink OL Input u_{dot}')
subplot(4,1,2); plot(time, u(:,1)); axis tight; grid on; 
ylabel('beta [rad]')
subplot(4,1,3); plot(time, uDot(:,2)); axis tight; grid on;
ylabel('TgDot_{ref} [Nm/s]'), 
subplot(4,1,4); plot(time, u(:,2)); axis tight; grid on;
ylabel('Tg_{ref} [Nm]'), xlabel('Time')

print('input_udot_u_OL', '-dpng');

% Table for state plots
xCellTable = {'omega_g', 'T_g', 'ThetaDot_s','omegaDot_r',...
'omegaDot_g','TDot_g'};
tableTestLPVxDot06 = array2table(xDot06,'VariableNames', xCellTable);
tableTestLPV = array2table(OutDataTestOL_LPV,'VariableNames', xCellTable);

% Output of linear model
figure; 
subplot(2,1,1); 
plot(time, tableTestLPVxDot06.omega_g,time,tableTestLPV.omega_g,'--');
axis tight; grid on;
ylabel('GenSpeed omegaG [rad/s]')
title('Output LPV model: Generator Speed \omega_g and Torque T_g')
subplot(2,1,2); 
plot(time, tableTestLPVxDot06.T_g, time, tableTestLPV.T_g,'--'); hold on
%plot(time, xk06(:,2))
legend('Discrete LPV x(t+\Delta t)','Discrete LPV x(t)')
axis tight; grid on;
ylabel('GenTorque T_G [Nm]')
print('output_y_OL', '-dpng');

% Output of linear and non linear model
figure; 
subplot(2,1,1); 
plot(time, tableTestLPVxDot06.omega_g, time, tableTestLPV.omega_g,'--',...
    time,xk06(:,1),'k-.')
axis tight; grid on;
ylabel('GenSpeed omegaG [rad/s]')
title('Output LPV model: Generator Speed \omega_g and Torque T_g')
subplot(2,1,2); 
plot(time, tableTestLPVxDot06.T_g, time, tableTestLPV.T_g,'--',time, xk06(:,2),'k-.')
legend('Discrete LPV x(t+\Delta t)','Discrete LPV x(t)', 'Nonlinear x(t)')
axis tight; grid on;
ylabel('GenTorque T_G [Nm]')
print('output_udotNonLin_OL', '-dpng');