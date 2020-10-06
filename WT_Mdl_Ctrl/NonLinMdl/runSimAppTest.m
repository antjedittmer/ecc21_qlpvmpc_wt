clear; clc; close all;  bdclose all; 
addpath(fullfile(fileparts(pwd),'datain'));
simMdlname = 'FAST_4TO24mpers_Bianchi'; %'FAST_4TO24mpers_SimApp'
%simMdelname = 'FAST_4TO24mpers_SimApp'; %
%simMdelname = 'FAST_4TO24mpers_Bianchi_TorqueCtrl';

load('workSpace.mat')
load('OutDataSweep.mat');
OutTable = array2table(OutData,'VariableNames', OutList);

idxTime = 1: height(OutTable);
OutTable = OutTable(idxTime,:);
OutTable.BldPitch1 = NaN(height(OutTable),1);
OutTable.BldPitch1(1:3:end) = pitch(:,1);

OutTable.Torque = NaN(height(OutTable),1);
OutTable.Torque(1:3:end) = Torque;

idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

param;
initModel5MWNREL;
initModel5MWNRELBianchi;

open_system(simMdlname);
sim(simMdlname);

varnames = {'Wind', 'RotSpeed', 'GenPwr', 'BldPitch1', 'GenTq', ...
    'NcIMUTAxs', 'NcIMUTAys'}; %%'Tower FA acc', 'Tower SS acc'};
OutTableTest = array2table(OutDataTest,'VariableNames',varnames);

%% Get actuator block to test influence
%actuatorBlock =  find_system(simMdelname,'Name','Actuators');
%Commented: 'through'/'off'


%% Prepare data for plots
% Calculate wind amplitude
cellChannelWind = {'Wind1Vel'};
OutList = OutTable.Properties.VariableNames;
vectAmpWindTest = sqrt(sum((vectWind.^2),2));

% Get color map for 'title legends'
testCaseStr = 'RotSpeed/TowerAcc';
cl = colormap('lines');
titleStr = ['Test Case ',testCaseStr,': {\color[rgb]{',num2str(cl(1,:)),'} FAST ',...
    '\color[rgb]{',num2str(cl(2,:)),'} Simulink} '];

%% Create output plots
figNo1 = 1; figNo2 = 2;

idxTime = 1: min(length(idxTime), height(OutTableTest));
time = OutTable.Time(idxTime);

defaultLineWidth = get(groot,'defaultLineLineWidth');
set(groot,'defaultLineLineWidth',1);

figure(figNo1);
axPlot(1) = subplot(3,1,1);
plot(time,vectAmpWind(idxTime),time,OutTableTest.Wind(idxTime),'--'); 
axis tight; grid on;
ylabel('Wind Amplitude (m/s)')
title(titleStr);

axPlot(2) = subplot(3,1,2);
plot(time,OutTable.RotSpeed(idxTime),time,OutTableTest.RotSpeed(idxTime),'--');
axis tight; grid on;
ylabel('RotSpeed (rpm)')

axPlot(3) = subplot(3,1,3);
plot(time,OutTable.GenPwr(idxTime)/1000,time,OutTableTest.GenPwr(idxTime)/1000,'--');
axis tight; grid on;
ylabel('GenPwr (MW)')
xlabel('Time (s)')
linkaxes(axPlot,'x');


figure(figNo2)
axPlot2(1) = subplot(4,1,1);
plot(time,OutTable.BldPitch1(idxTime),'.', time,OutTableTest.BldPitch1(idxTime),'--'); 
axis tight; grid on;
ylabel('BldPitch1 (rad)')
title(titleStr);

axPlot2(2) = subplot(4,1,2);
plot(time,OutTable.GenTq(idxTime),time,OutTableTest.GenTq(idxTime)/10^3);
axis tight; grid on;
ylabel('GenTq (Nm)') 

axPlot2(3) = subplot(4,1,3);
plot(time,OutTable.NcIMUTAxs(idxTime),time,OutTableTest.NcIMUTAxs(idxTime),'--');
axis tight; grid on;
ylabel('Tower FA (m/s)')

axPlot2(4) = subplot(4,1,4);
plot(time,OutTable.NcIMUTAys(idxTime),time,OutTableTest.NcIMUTAys(idxTime),'--');
axis tight; grid on;
ylabel('Tower SS (m/s)')
xlabel('Time (s)')
linkaxes(axPlot2,'x');
set(groot,'defaultLineLineWidth',defaultLineWidth);

%% Get incremental torque characteristic TheadG = B
vectAmpWind(idxTime);
lengthWindSteps = 5:24; 

vecIdx = (2500:2500:length(idxTime))- 50;
genSpeed_rad_s = OutTable.GenSpeed(vecIdx) *pi/30; 
genTq_kNm = OutTable.GenTq(vecIdx);

ratioGenTqSpeed = genTq_kNm ./ genSpeed_rad_s;
figure; plot(genTq_kNm ./ genSpeed_rad_s);
grid on; axis tight;
title('genTq_{kNm}/genSpeed_{rad,s}')

%[ratioGenTqSpeed(1:20),(5:24)']
Be_5To14 = [0.0664, 0.1109,0.1598,0.2086,0.2358,0.2623,0.3079,0.3574,0.3576,0.3577];


%% plot measured x = [y,xdot]
if ~exist('xk','var')==1
    return
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