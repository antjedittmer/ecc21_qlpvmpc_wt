%% get values
clear;
clear ABgL_WT;
close all;

% set paths
dataInPath = fullfile(fileparts(pwd),'dataIn');
addpath(dataInPath);
addpath(fullfile(pwd,'dataOut'));
addpath(genpath('C:\Projects\FASTTool'));
addpath('subfunctionsMPC')
addpath('MDlWECS')
onlyInit = 1;
dataInPath = fullfile(fileparts(pwd),'dataIn');
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable'); %#ok<NASGU>
DT = 0.008;

% if onlyInit == 1
%     return;
% end

%initModel5MWNREL;
initModel5MWNREL;

% set local fst file
FAST_InputFileName = fullfile(pwd,'inputfiles', 'FAST.fst');
selectSteadyState = 1; % 0: OpenLoop FAST, 1: PI CL FAST, 2:PI CL Simulink

load('workSpace.mat')
load('OutDataSweep.mat','OutData');
%load('OutDataSweep.mat','OutData');
%load(fullfile(dataInPath,'OutDataWind18NTM.mat'),'OutData');
%load('OutDataWind11To14NTM.mat','OutData');

OutTable = array2table(OutData,'VariableNames', OutList);
idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));
OutTable.Wind1VelX = vectAmpWind;

%% LOok up table for stable states
% LUT_V = Lin_points.V; %4:25;
% LUT_omega = Lin_points.RSpeed; %[0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)]; %Lin_points.RSpeed
% LUT_beta =  Lin_points.Pitch; % [zeros(1,8), 0.0642    0.1126    0.1493    0.1808    0.2093    0.2356    0.2602...    0.2833    0.3054    0.3269    0.3477    0.3676    0.3868    0.4056];
% LUT_torque = Lin_points.Torque; [5.0870    7.9450   11.4300   15.5300   20.2400   25.5200   31.3500  ,repmat(43.0900,1,15)];

if selectSteadyState == 2 % select CL PI Simulink steady state
    load('OutTableTest1')
    %OutTableTest1 = OutTableTest2;
    %clear('OutTableTest1');
    %OutTableTest1 = OutTable;
    lenTable = min(height(OutTableTest1), height(OutTable));
    
    OutTable = OutTable(1:lenTable,:);
    vectAmpWind = vectAmpWind(1:lenTable);
    
    OutTableTest1 = OutTableTest1(1:lenTable,:);
    OutTableTest1.Time = OutTable.Time;
    OutTableTest1.GenTq = OutTable.GenTq; %N ->k Nm (same unit as Linpoints)
    OutTableTest1.Wind1VelX  = vectAmpWind;
    
    OutTable = OutTableTest1;
    vectAmpWind = OutTable.Wind1VelX;

end


if selectSteadyState > 0 % overwrite default FAST OL data   
    load(fullfile(dataInPath,'Lin_points_PICL.mat'),'Lin_points');
end
 

%% Physical parameters
% Out: omegaG, Tg
% States: ThetaS, omegaR, omegaG, Tg
% In: beta_ref, Tg_ref;
% Scheduling param rho: V, beta_ref, omegaG/Ng (as approximation of omegaR)

%vectAmpWind = OutTable.Wind1VelX(20/DT:end);
Vbar = OutTable.Wind1VelX(1);
rb = 63; % m blade radius
G = 97;
Ng = G;
tau_g = 0.01;
omegabar = interp1(Lin_points.V,Lin_points.RSpeed,Vbar,'spline');
Lambda_bar = omegabar*wecs.rb/Vbar;
beta_bar = interp1(Lin_points.V,Lin_points.Pitch,Vbar,'spline');
schedParam = [omegabar, beta_bar, Vbar];

r(1) = interp1(Lin_points.V,Lin_points.RSpeed,Vbar,'spline') * Ng; %
r(2) = interp1(Lin_points.V,Lin_points.Torque * 10^3,Vbar,'spline');

xk = [r(1), r(2), zeros(1,4),0]';

r7(1) = r(2);
r7(2) = interp1(Lin_points.V,Lin_points.Pitch,Vbar,'spline'); %
xk07 =  [r7(1), r7(2), zeros(1,5),0]';

nu = 2;

N = 40;
Xko = kron(ones(N,1),xk);
Xko(7:7:end)= 0;

Xk07 = kron(ones(N,1),xk07);
Xk07(8:8:end)= 0;

mu_old = zeros(2*N*nu,1);
Rotor_cQ = Rotor_cQ(:,1:size(Rotor_Lambda,1)); lambdaVec = Rotor_Lambda;
pitchVec = Rotor_Pitch * pi/180;

U = zeros(2*N,1);
cQ = Rotor_cQ;

idxWind_11to14 = vectAmpWind >= 10 &  15 > vectAmpWind;
vectAmpWind_11to14 = vectAmpWind(idxWind_11to14);
time_11to14 = 0:DT: DT*length(vectAmpWind_11to14)-DT;

if onlyInit == 1
    return
end


% Test function
[u,Xko,Unew,mu,us,tictoc]= qLPVMPC_WECS_matlab(r,xk,Xko,mu_old,lambdaVec,pitchVec,cQ,schedParam,U,DT);


fprintf('Reference value from steady state r: '); fprintf('%2.3f ',r); fprintf('\n');
fprintf('Calculated input beta [rad] %2.3f\n',u(1));
fprintf('Calculated input T_g [Nm] %2.3f\n',u(2));
fprintf('Estimated states: '); fprintf('%2.3f ', Xko(1:6)); fprintf('\n');
fprintf('Estimated states next step: '); fprintf('%2.3f ', Xko((1:6) +6)); fprintf('\n');

%% Run Matlab simulation in 10 steps
% Call the Matlab simulation 10 times with constant scheduling parameters
load('OutDataTestCL_LPV_01','xk06','UVecSim','XVecSim','u_uss','tictoc_LPVMPC','q_','r_','p');
% Initialize vectors and matrices for simulation
Nsim = 10;
U = zeros(2*N,1);
xDdot_Ve = NaN(Nsim+1,length(xk));
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
xDdot_Ve = NaN(Nsim+1,7);
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
load('OutDataTestCL_LPV_07', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B','UVecSim','XVecSim','u_uss','tictoc_LPVMPC','q_','r_','p');
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
u1Cell  = {'beta[deg]', 'T_g [Nm]'};

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
load('OutDataTestOL_LPV_windChange_200sec', 'A', 'B', 'u', 'uDot', 'wind', 'xDot06', 'xk06', 'xkLPV');% 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B');

% Plot input
time = OutTable.Time(1: length(uDot));
figure; 
subplot(4,1,4); plot(time, uDot(:,2)); axis tight; grid on; 
ylabel('\Delta \beta [rad]'); xlabel('time (s)')
title('Simulink OL Input u_{dot}')
subplot(4,1,3); plot(time, u(:,2)); axis tight; grid on; 
ylabel('pitch \beta [rad]')
subplot(4,1,2); plot(time, uDot(:,3)); axis tight; grid on;
ylabel('\Delta Tg_{ref} [Nm]'), 
subplot(4,1,1); plot(time, u(:,3)); axis tight; grid on;
ylabel('torque Tg_{ref} [Nm]'), 
print('input_udot_u_OL', '-dpng');

P_g = xk06(:,1) .* xk06(:,2)/10^6;

figure; 
subplot(4,1,1); plot(time, wind); axis tight; grid on; 
ylabel('wind [m/s]'); %xlabel('time (s)')
title('Windsweep Variable-speed Variable-pitch Wind Turbine');
subplot(4,1,2); plot(time,u(:,3)); axis tight; grid on;
ylabel('torque Tg_{ref} [Nm]'); 
subplot(4,1,3); plot(time, u(:,2)*180/pi); axis tight; grid on; 
ylabel('pitch \beta [deg]')
subplot(4,1,4); plot(time, P_g); axis tight; grid on; 
ylabel('power P_g [MW]')
xlabel('time [s]');

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.25]);

set(findall(gcf,'-property','FontSize'),'FontSize',11)
print('windsweep_u_y_OL', '-dpng');

% Table for state plots
idxTime =  1 : min([length(time), length(xkLPV), length(xDot06)]);
xCellTable = {'omega_g', 'T_g', 'ThetaDot_s','omegaDot_r',...
'omegaDot_g','TDot_g'};
tableTestLPVxDot06 = array2table(xDot06,'VariableNames', xCellTable);
tableTestLPV = array2table(xkLPV,'VariableNames', xCellTable);
tableTestLPVxDot06 = tableTestLPVxDot06(idxTime,:);
tableTestLPV = tableTestLPV(idxTime,:);

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
