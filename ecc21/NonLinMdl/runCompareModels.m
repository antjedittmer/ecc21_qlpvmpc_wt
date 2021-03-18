function runCompareModels(strWindType,loadData,yAxCell)
% runCompareModels compares two Simulink models (Mdl1(5 states) and Mdl2 (8 
% states) with 

if ~nargin || isempty(strWindType)
    strWindType = 'Sweep'; 
end

if nargin < 2 || isempty(loadData)
   loadData = 1; %load simulation output data if available;
end

if nargin < 3 || isempty(yAxCell)
    yAxCell = {'wind V [m/s]', 'GenTq T_g [kNm]', 'Pitch \beta [°]', 'RotSpd \omega_r [rpm]',...
    'GenPwr P_g [MW]','Twr_{FA} y_t [m/s^2]'};
end

%% Initialize workspace and path
%lear; clc; close all;  
% bdclose all; 

% Set path to inputdata directory
mainDir = fullfile(pwd);
dataInDir = fullfile(fileparts(pwd),'dataIn');
addpath(dataInDir);

% Set path to figure directory
%figDir = fileparts(mainDir); %'html','figDir');
figDir = fullfile(fileparts(mainDir),'figDir');

if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory
dataDirOut = 'dataOut';
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Provide names of FAST simulation data to be loaded
%'NTW18';%'Sweep'; %normal turbulance

% Provide names of Simulink models to be run
simMdlname1 = 'test_SimulinkMdl1_Baseline'; %
%simMdlname2 = 'test_SimulinkMdl3LPV_Baseline'; 
simMdlname2 = 'test_SimulinkMdl2_Baseline'; 

if strcmp(strWindType,'Sweep') == 1
    % sweep from 4 to 25 in steps
    %workspaceSimulationMat = 'workSpace.mat';
    outDataSimulationMat = 'OutDataSweep.mat';
    strFig = '';
    testCaseStr = 'Wind Sweep';
else
    % wind with average 18 m/s
    %workspaceSimulationMat = 'workSpace18NTM.mat';
    outDataSimulationMat = 'OutDataWind18NTM.mat';
    strFig = 'NTW18';
    testCaseStr = 'Wind, mean 18 m/s';
end

%% Load data from FAST run
%load(fullfile(dataInDir,workspaceSimulationMat)) %#ok<LOAD> ADi ToDo: Clean this up later: Only load variables needed for simulation
% DT = 0.008;
load(outDataSimulationMat ,'OutTable');%#ok<NASGU>

% LUT tables loaded have these values
% LUT_torque2 = [5.0066    8.8294   13.5530   18.9143   23.8603   29.5200   36.9452,...
%    repmat(43.5290,1,15)];
% LUT_omega2 = [0.7776    0.8209    0.8742    0.9348    1.0430    1.1600  ...
%     1.2369    1.2555    1.2548    1.2546    1.2546    1.2545    1.2546    1.2545    1.2545    1.2545    1.2545 ...
%     1.2545    1.2545    1.2545 1.2545    1.2545];

% Assign data from FAST to variables
% OutTable = array2table(OutData,'VariableNames', OutList);
idxTime = 1: height(OutTable);
OutTable = OutTable(idxTime,:);
% OutTable.BlPitch1 = NaN(height(OutTable),1);
% OutTable.BldPitch1(1:3:end) = pitch(:,1);  %#ok<NODEF> Loaded in Matfile 
OutTable.BldPitch1 = OutTable.BlPitch1;

% OutTable.Torque = NaN(height(OutTable),1);
% OutTable.Torque(1:3:end) = Torque;
OutTable.Torque = OutTable.GenTq;

idxWind = contains(OutTable.Properties.VariableNames, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

% assign data for LPV model
%initModel5MWNREL;

% idxWind_11to14 = vectAmpWind > 10 &  14 > vectAmpWind;
% vectAmpWind_11to14 = vectAmpWind(idxWind_11to14);
% time_11to14 = 0:DT: DT*length(vectAmpWind_11to14)-DT;

% Initialize constants for Simulink model
% param;
%initModel5MWNREL;

% Lin_points.RSpeed = LUT_omega2;
% Lin_points.Torque = LUT_torque2;
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', ...
    'NcIMUTAxs', 'NcIMUTAys'}; %%'Tower FA acc', 'Tower SS acc'};

%% Run simulations or load data
% Run simulation for two models or load mat files if available
matFileOutTableTest1 = fullfile(dataDirOut,['OutTableTest1',strFig,'.mat']);
if exist(matFileOutTableTest1,'file') == 2 && loadData
    load(matFileOutTableTest1,'OutTableTest1');
else
    WS = 'base';
    % assignin(WS,'DT',DT)
    assignin(WS,'OutTable',OutTable)
    open_system(simMdlname1);
    %Ng = wecs.Ng;
    out1 = sim(simMdlname1,'SrcWorkspace',WS);
    OutTableTest1 = array2table(out1.OutDataTest,'VariableNames',varnames);
    save(matFileOutTableTest1,'OutTableTest1');
end

matFileOutTableTest2 = fullfile(dataDirOut,['OutTableTest2',strFig,'.mat']);
if exist(matFileOutTableTest2,'file') == 2 && loadData
    load(matFileOutTableTest2,'OutTableTest2');
else
    open_system(simMdlname2);
    out2 = sim(simMdlname2,'SrcWorkspace',WS);
    if size(out2.OutDataTest,2) == 9 % Pitch commands of all actuators
        out2OutDataTest = out2.OutDataTest(:,[1:4,7:9]);
    else
         out2OutDataTest = out2.OutDataTest;
    end
    OutTableTest2 = array2table(out2OutDataTest,'VariableNames',varnames);
    save(matFileOutTableTest2,'OutTableTest2');
end

% Get actuator block to test influence
% actuatorBlock =  find_system(simMdelname,'Name','Actuators');
% Commented: 'through'/'off'

%% Create output plots
%Prepare data for plots
% Calculate wind amplitude
%cellChannelWind = {'Wind1Vel'};
%OutList = OutTable.Properties.VariableNames;
%vectAmpWindTest = sqrt(sum((vectWind.^2),2));

% Get color map for 'title legends'
cl = colormap('lines');
titleStr = [testCaseStr,': Mdl1: Rot+Twr {\color[rgb]{',num2str(cl(1,:)),'}Mdl2: Rot,Twr,Bld+Act ',...
    '\color[rgb]{',num2str(cl(2,:)),'}FAST} '];

%Create output plots
figNo1 = 1; figNo2 = 2;

idxTime = 1: min([length(idxTime), height(OutTableTest1),height(OutTableTest2)]);
time = OutTable.Time(idxTime);

idxTime = time > 15;
time = time(idxTime);

% defaultLineWidth = get(groot,'defaultLineLineWidth');
% set(groot,'defaultLineLineWidth',1);
% set(groot,'defaultLineLineWidth',defaultLineWidth );

figure(figNo1);
axPlot(1) = subplot(3,1,1);
plot(time,vectAmpWind(idxTime),time,OutTableTest1.Wind(idxTime),'--'); 
axis tight; grid on;
ylabel('wind [m/s]')
title(titleStr);

axPlot(2) = subplot(3,1,2);
plot(time,OutTableTest2.RotSpeed(idxTime), time,OutTable.RotSpeed(idxTime),time,OutTableTest1.RotSpeed(idxTime),'k--');
axis tight; grid on;
ylabel('RotSpd \omega_r [rpm]')

axPlot(3) = subplot(3,1,3);
plot(time,OutTableTest2.GenPwr(idxTime)/1000,time,OutTable.GenPwr(idxTime)/1000,time,OutTableTest1.GenPwr(idxTime)/1000,'k--');
axis tight; grid on;
ylabel('GenPwr P_g [MW]')
xlabel('time [s]')
linkaxes(axPlot,'x');

set(findall(gcf,'-property','FontSize'),'FontSize',11)
print(fullfile(figDir,['cmpTimeDomain_Wind',strFig]), '-dpng');

figure(figNo2)
axPlot2(1) = subplot(3,1,1);
plot(time,OutTableTest2.BlPitch1(idxTime),time,OutTable.BlPitch1(idxTime),time,OutTableTest1.BlPitch1(idxTime),'k--'); 
axis tight; grid on;
ylabel('BldPitch1 \beta [rad]')
title(titleStr);

axPlot2(2) = subplot(3,1,2);
plot(time,OutTableTest2.GenTq(idxTime)/10^3,time,OutTable.GenTq(idxTime),time,OutTableTest1.GenTq(idxTime)/10^3,'k--');
axis tight; grid on;
ylabel('GenTq T_g [kNm]') 

axPlot2(3) = subplot(3,1,3);
plot(time,OutTableTest2.NcIMUTAxs(idxTime),time,OutTable.NcIMUTAxs(idxTime),time,OutTableTest1.NcIMUTAxs(idxTime),'k--');
axis tight; grid on;
ylabel('Twr_{FA} y_t[m/s^2]')
xlabel('time [s]')
set(findall(gcf,'-property','FontSize'),'FontSize',11)
print(fullfile(figDir,['cmpTimeDomain_BlPitchTqTwr',strFig]), '-dpng');
linkaxes(axPlot2,'x');


%% Figure for paper
figure(figNo2+1)
axPlotAll(1) = subplot(6,1,1);
plot(time,vectAmpWind(idxTime),time,OutTableTest1.Wind(idxTime),'--'); 
axis tight; grid on;
ylabel(yAxCell{1}); %'wind [m/s]')
title(titleStr);
  
axPlotAll(4) = subplot(6,1,4);
plot(time,OutTableTest2.RotSpeed(idxTime), time,OutTable.RotSpeed(idxTime),time,OutTableTest1.RotSpeed(idxTime),'k--');
axis tight; grid on;
ylabel(yAxCell{4}); %'RotSpd \omega_r [rpm]')

axPlotAll(2) = subplot(6,1,2);
plot(time,OutTableTest2.GenTq(idxTime)/10^3,time,OutTable.GenTq(idxTime),time,OutTableTest1.GenTq(idxTime)/10^3,'k--');
axis tight;
posAxis = axis;
axis([posAxis(1:2), min(43,posAxis(3)), 44]);
ylabel(yAxCell{3}); %'T_g [kNm]')
grid on; % axis tight;
%grid on;
ylabel(yAxCell{2}); %'GenTq T_g [Nm]') 

axPlotAll(3) = subplot(6,1,3);
plot(time,OutTableTest2.BlPitch1(idxTime),time,OutTable.BlPitch1(idxTime),time,OutTableTest1.BlPitch1(idxTime),'k--'); 
axis tight; grid on;
ylabel(yAxCell{3}); %'Pitch \beta [°]')

axPlotAll(5) = subplot(6,1,5);
plot(time,OutTableTest2.NcIMUTAxs(idxTime),time,OutTable.NcIMUTAxs(idxTime),time,OutTableTest1.NcIMUTAxs(idxTime),'k--');
axis tight; grid on;
ylabel(yAxCell{6}); %'Twr_{FA} y_t[m/s^2]')

axPlotAll(6) = subplot(6,1,6);
plot(time,OutTableTest2.GenPwr(idxTime)/1000,time,OutTable.GenPwr(idxTime)/1000,time,OutTableTest1.GenPwr(idxTime)/1000,'k--');
axis tight; grid on;
ylabel(yAxCell{5}); %'GenPwr P_g [MW]')

xlabel('time [s]')

linkaxes(axPlotAll,'x');

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.5]);
print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-dpng');
print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');

% axPlot2(4) = subplot(4,1,4);
% plot(time,OutTableTest2.NcIMUTAys(idxTime),time,OutTable.NcIMUTAys(idxTime),time,OutTableTest1.NcIMUTAys(idxTime),'k--');
% axis tight; grid on;
% ylabel('Twr_{SW} x_t (m/s^2)')
% xlabel('time (s)')
% linkaxes(axPlot2,'x');
% set(groot,'defaultLineLineWidth',defaultLineWidth);

return;

%% Get incremental torque characteristic TheadG = B
lengthWindSteps = 5:24;  %#ok<UNRCH> This was used for debugging - ADi ToDo delete or move to another file

vecIdx = (2500:2500:length(idxTime))- 50;
%vecIdx = (2500:2500: height(OutTable)) -100

wind_Vec = round(vectAmpWind(vecIdx));

genSpeed_rad_s = OutTable.GenSpeed(vecIdx) * pi/30; 
genTq_kNm = OutTable.GenTq(vecIdx);
genPwr_kNm = OutTable.GenPwr(vecIdx);

genSpeed_rad_s_Mdl1 = OutTableTest1.RotSpeed(vecIdx) *pi/30 * 97; 
genTq_kNm_Mdl1 = OutTableTest1.GenTq(vecIdx)/1000;

genSpeed_rad_s_Mdl2 = OutTableTest2.RotSpeed(vecIdx) *pi/30 * 97; 
genTq_kNm_Mdl2 = OutTableTest2.GenTq(vecIdx)/1000;

try
LUT_Torque = interp1(Lin.V, Lin.Torque, wind_Vec);
catch
    return
end
LUT_GenSpeed = interp1(Lin.V, Lin.RSpeed, wind_Vec) * 97;

figure(figNo2+1)
subplot(3,1,1)
plot(wind_Vec,LUT_GenSpeed,'o', wind_Vec,genSpeed_rad_s,'*', wind_Vec,genSpeed_rad_s_Mdl2,'kx')
title('Reference Values: Steady-state Generator Speed and Torque')
ylabel('GenSpeed \omega_g [rad/s]')
axis tight; grid on;

subplot(3,1,2)
plot(wind_Vec,LUT_Torque,'o', wind_Vec,genTq_kNm,'*', wind_Vec,genTq_kNm_Mdl2,'kx')
ylabel('GenTq T_g [kNm]')
axis tight; grid on;

subplot(3,1,3)
plot(wind_Vec,LUT_Torque .* LUT_GenSpeed,'o', ...
    wind_Vec,genTq_kNm .* genSpeed_rad_s,'*', ...
    wind_Vec,genSpeed_rad_s_Mdl2 .* genTq_kNm_Mdl2,'kx')
ylabel('GenPwr P_g [kW]')
xlabel('Wind speed V [m/s]')
axis tight; grid on;
legend('FAST LUT','FAST Sim','Simulink MDL2 Sim', 'Location','SouthEast')
print(fullfile(figDir,['cmpReferences_LUT_FAST_Simulink',strFig]), '-dpng');

%ratioGenTqSpeed = genTq_kNm ./ genSpeed_rad_s;

figure(figNo2+2)
subplot(2,2,1); plot(OutTable.GenSpeed,OutTable.GenTq);
ylabel('GenTq T_{g,ref} [kNm]');
grid on; axis tight;
subplot(2,2,3); plot(OutTable.GenSpeed,OutTable.BlPitch1);
ylabel('Pitch \beta_ref [kNm]');
xlabel('GenSpeed \omega_g [rpm]')
grid on; axis tight;
subplot(2,2,2); plot(vectAmpWind,OutTable.GenTq);
grid on; axis tight;
subplot(2,2,4); plot(vectAmpWind,OutTable.BlPitch1);
grid on; axis tight;
xlabel('Wind Speed V [m/s]')

%title('genTq_{kNm}/genSpeed_{rad,s}')

%[ratioGenTqSpeed(1:20),(5:24)']
Be_5To14 = [0.0664, 0.1109,0.1598,0.2086,0.2358,0.2623,0.3079,0.3574,0.3576,0.3577];


%% plot measured x = [y,xdot]
if ~exist('xk','var')==1
    return
end

if min(size(xk)) == 1
    return;
end

xest = xk;
xCell = {'ThetaDot_s [rad]','omegaDot_r [rad]',...
    'omegaDot_g [rad]','TDot_g','BetaDot_r','V_dot'}; %'V_dot'

figure(4); 
subplot(3,1,1); plot(xest(:,1))
grid on; axis tight
legend('omega_g [rad/s]')
subplot(3,1,2);
plot(xest(:,3)); hold on; plot(xest(:,4),'k')
grid on; axis tight
legend('beta [deg]','V [m/s]');
subplot(3,1,3); plot(xest(:,2));legend('T_g [Nm]'); 
grid on; axis tight

figure(5);
for idx = 1:3
    subplot(3,1,idx);
    plot(xest(:,idx+4))
    axis tight; grid on;
    ylabel(xCell{idx});
end

figure(6);
for idx = 1:3
    subplot(3,1,idx);
    plot(xest(:,idx+4+3))
    axis tight; grid on;
    ylabel(xCell{idx+3})
end