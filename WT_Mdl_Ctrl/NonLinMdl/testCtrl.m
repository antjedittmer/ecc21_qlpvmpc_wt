%% get values
clear;
close all;
addpath(fullfile(fileparts(pwd),'dataIn'));
initModel5MWNRELBianchi;

load('workSpace.mat')
load('OutDataSweep.mat');
OutTable = array2table(OutData,'VariableNames', OutList);
idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

%% LOok up table for stable states
LUT_V = 4:25;
LUT_omega = [0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)];
LUT_beta = [zeros(1,8),0.0642    0.1126    0.1493    0.1808    0.2093    0.2356    0.2602...
    0.2833    0.3054    0.3269    0.3477    0.3676    0.3868    0.4056];
LUT_torque = [5.0870    7.9450   11.4300   15.5300   20.2400   25.5200   31.3500  ,...
    repmat(43.0900,1,15)];

%% Physical parameters
% Out: omegaG, Tg
% States: ThetaS, omegaR, omegaG, Tg
% In: beta_ref, Tg_ref;
% Scheduling param rho: V, beta_ref, omegaG/Ng (as approximation of omegaR)
Vbar = 4;
rb = 63; % m blade radius
G = 97;
omegabar = interp1(LUT_V,LUT_omega,Vbar);
Lambda_bar = omegabar*rb/Vbar;
beta_bar = interp1(LUT_V,LUT_beta,Vbar)/pi*180;
schedParam = [omegabar, beta_bar/180* pi, Vbar];

r(1) = interp1(LUT_V,LUT_omega,Vbar) * Ng; %
r(2) = interp1(LUT_V,LUT_torque,Vbar) * 1000;

xk = [r(1), r(2), zeros(1,4)]';
nu = 2;

N = 40;
Xko = kron(ones(N,1),xk);
mu_old = zeros(2*N*nu,1);
Rotor_cQ = Rotor_cQ(:,1:size(Rotor_Lambda,1)); lambdaVec = Rotor_Lambda;
pitchVec = Rotor_Pitch;

U = zeros(2*N,1);
cQ = Rotor_cQ;

% Test function
[u,Xko,Unew,mu,us,tictoc]= qLPVMPC_WECS_matlab(r,xk,Xko,mu_old,lambdaVec,pitchVec,cQ,schedParam,U,DT);
fprintf('Reference value from steady state r: '); fprintf('%2.3f ',r); fprintf('\n');
fprintf('Calculated input beta [rad] %2.3f\n',u(1));
fprintf('Calculated input T_g [Nm] %2.3f\n',u(2));
fprintf('Estimated states: '); fprintf('%2.3f ', Xko(1:6)); fprintf('\n');
fprintf('Estimated states next step: '); fprintf('%2.3f ', Xko((1:6) +6)); fprintf('\n');

%% Run Matlab simulation in 10 steps
% Call the Matlab simulation 10 times with constant scheduling parameters

% Initialize vectors and matrices for simulation
Nsim = 10;
U = zeros(2*N,1);
xDdot_Ve = NaN(Nsim+1,6);
uvec = NaN(Nsim+1,2);
Uvec = NaN(Nsim+1,80);

% Values at time 0 for simulation
xDdot_Ve(1,:) = xk;
Uvec(1,:) = U;
uvec(1,:) = [schedParam(2), Xko(2)];
nxa = length(xk);

for i = 1:Nsim
    [u,Xko,U,mu,us,tictoc]= qLPVMPC_WECS_matlab(r,xk,Xko,mu_old,lambdaVec,pitchVec,cQ,schedParam,U,DT);
    xk = Xko(nxa+1:2*nxa);
    xDdot_Ve(i+1,:) = xk;
    Uvec(i+1,:) = U;
    uvec(i+1,:) = u;
end
fprintf('Calculated input beta [rad] %2.3f\n',u(1));
fprintf('Calculated input T_g [Nm] %2.3f\n',u(2));

time = 0: DT:(Nsim -1)*DT;
titleStr1 = 'Input Udot Matlab Simulation';
titleStr2 = ['State Estimates Matlab Simulation, const speed ',num2str(Vbar)];
nrPlot1 = 1;
plotInputStateTrajectories(Uvec,xDdot_Ve,time, titleStr1, titleStr2,nrPlot1);

%% Plot one time step simulation results:
% load('OutDataTestCL_LPV_01', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B');
load('testVar','r','xk','Xk','mu_old','lambdaVec','pitchVec','cQ','schedParam','DT');
fprintf('Reference value from Simulink r: '); fprintf('%2.3f ',r); fprintf('\n');

% Initialize vectors and matrices for simulation
Nsim = 10;
U = zeros(2*N,1);
xDdot_Ve = NaN(Nsim+1,6);
uvec = NaN(Nsim+1,2);
Uvec = NaN(Nsim+1,80);

% Values at time 0 for simulation
xDdot_Ve(1,:) = xk;
Uvec(1,:) = U;
uvec(1,:) = [schedParam(2), Xk(2)];

for i = 1:Nsim
    [u,Xk,U,mu,us,tictoc]= qLPVMPC_WECS_matlab(r,xk,Xk,mu_old,lambdaVec,pitchVec,cQ,schedParam,U,DT);
    xk = Xk(nxa+1:2*nxa);
    xDdot_Ve(i+1,:) = xk;
    Uvec(i+1,:) = U;
    uvec(i+1,:) = u;
end
fprintf('Calculated input beta [rad] %2.3f\n',u(1));
fprintf('Calculated input T_g [Nm] %2.3f\n',u(2));

time = 0: DT:(Nsim -1)*DT;
titleStr1 = 'Input Udot Matlab Simulation';
titleStr2 = ['State Estimates Matlab Simulation, const speed ',num2str(Vbar)];
nrPlot1 = 3;
plotInputStateTrajectories(Uvec,xDdot_Ve,time, titleStr1, titleStr2,nrPlot1);

%% Same plots for Simulink simulations
load('OutDataTestCL_LPV_06', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B','UVecSim','XVecSim','u_uss','tictoc_LPVMPC','q_','r_','p');
nrPlot1 = 5;
Uvec =  UVecSim; 
xDdot_Ve = xk06;
Nsim = length(u_uss);
time = OutTable.Time(1: Nsim);
plotInputStateTrajectories(Uvec,xDdot_Ve,time, titleStr1, titleStr2,nrPlot1);

%% Plot results from closed loop simulation
%load('OutDataTestCL_LPV_01', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B','UVecSim','XVecSim','u_uss','tictoc_LPVMPC','q_','r_','p');
u_uss_plot = squeeze(u_uss);
u_uss_plot = diag([180/pi,1,180/pi,1]) *  u_uss_plot';
u1Cell  = {'beta[deg]','T_g [Nm]'};

wind = vectAmpWind;
timeEnd = Nsim;
Nsim = length(u_uss);
time = OutTable.Time(1: length(u_uss));

qVec = deblank(sprintf('%d ',q_));
rVec = deblank(sprintf('%d ', r_));
pVec = deblank(sprintf('%d', p));

figure(11);
subplot(3,1,1);
plot(time(1:timeEnd),  wind(1:timeEnd))
ylabel('wind [m/s]')

title({'Simulation run with weights', ...
    ['q_ = [', qVec,'], r = [',rVec,'], p = ',pVec,'* q']})
R_ = kron(eye(N),diag([1,100]));
axis tight; grid on;
for idx = 1:2
    
    subplot(3,1,idx+1);
    plot(time(1:timeEnd), u_uss_plot(idx,1:timeEnd)); hold on;
    plot(time(1:timeEnd), u_uss_plot(idx+2,1:timeEnd));
    axis tight; grid on;
    ylabel(u1Cell{idx});
end
legend('u_{opt}','u');
hold off


%% Load data from OL simulation
load('OutDataTestOL_LPV', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B');

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
% tableTestLPVxDot06 = array2table(xDot06,'VariableNames', xCellTable);
% tableTestLPV = array2table(OutDataTestOL_LPV,'VariableNames', xCellTable);


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
