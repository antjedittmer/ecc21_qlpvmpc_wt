clear; clc; close all;  bdclose all; 

addpath(fullfile(fileparts(pwd),'dataIn'));
simMdlname1 = 'FAST_4TO24mpers_SimApp'; %
simMdlname2 = 'FAST_4TO24mpers_Bianchi'; %'FAST_4TO24mpers_SimApp'
% simMdelname = 'FAST_4TO24mpers_Bianchi_TorqueCtrl';

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

idxWind_11to14 = vectAmpWind > 10 &  14 > vectAmpWind;
vectAmpWind_11to14 = vectAmpWind(idxWind_11to14);
time_11to14 = 0:DT: DT*length(vectAmpWind_11to14)-DT;


param;
initModel5MWNREL;
initModel5MWNRELBianchi;
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'BlPitch1', 'GenTq', ...
    'NcIMUTAxs', 'NcIMUTAys'}; %%'Tower FA acc', 'Tower SS acc'};

open_system(simMdlname1);
sim(simMdlname1);
OutTableTest1 = array2table(OutDataTest,'VariableNames',varnames);

open_system(simMdlname2);
sim(simMdlname2);
OutTableTest2 = array2table(OutDataTest,'VariableNames',varnames);

%% Get actuator block to test influence
%actuatorBlock =  find_system(simMdelname,'Name','Actuators');
%Commented: 'through'/'off'


%% Prepare data for plots
% Calculate wind amplitude
cellChannelWind = {'Wind1Vel'};
OutList = OutTable.Properties.VariableNames;
vectAmpWindTest = sqrt(sum((vectWind.^2),2));

% Get color map for 'title legends'
testCaseStr = 'Wind Sweep';
cl = colormap('lines');
titleStr = [testCaseStr,': Mdl1: Rot+Twr {\color[rgb]{',num2str(cl(1,:)),'}Mdl2: Rot,Twr,Bld+Act ',...
    '\color[rgb]{',num2str(cl(2,:)),'}FAST} '];


%% Create output plots
figNo1 = 1; figNo2 = 2;

idxTime = 1: min([length(idxTime), height(OutTableTest1),height(OutTableTest2)]);
time = OutTable.Time(idxTime);

defaultLineWidth = get(groot,'defaultLineLineWidth');
set(groot,'defaultLineLineWidth',1);

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
print('cmpTimeDomain_Wind', '-dpng');


figure(figNo2)
axPlot2(1) = subplot(3,1,1);
plot(time,OutTableTest2.BlPitch1(idxTime),time,OutTable.BlPitch1(idxTime),time,OutTableTest1.BlPitch1(idxTime),'k--'); 
axis tight; grid on;
ylabel('BldPitch1 \beta [rad]')
title(titleStr);

axPlot2(2) = subplot(3,1,2);
plot(time,OutTableTest2.GenTq(idxTime)/10^3,time,OutTable.GenTq(idxTime),time,OutTableTest1.GenTq(idxTime)/10^3,'k--');
axis tight; grid on;
ylabel('GenTq \omega_g [Nm]') 

axPlot2(3) = subplot(3,1,3);
plot(time,OutTableTest2.NcIMUTAxs(idxTime),time,OutTable.NcIMUTAxs(idxTime),time,OutTableTest1.NcIMUTAxs(idxTime),'k--');
axis tight; grid on;
ylabel('Twr_{FA} y_t[m/s^2]')
xlabel('time [s]')
set(findall(gcf,'-property','FontSize'),'FontSize',11)
print('cmpTimeDomain_BlPitchTqTwr', '-dpng');

% axPlot2(4) = subplot(4,1,4);
% plot(time,OutTableTest2.NcIMUTAys(idxTime),time,OutTable.NcIMUTAys(idxTime),time,OutTableTest1.NcIMUTAys(idxTime),'k--');
% axis tight; grid on;
% ylabel('Twr_{SW} x_t (m/s^2)')
% xlabel('time (s)')
% linkaxes(axPlot2,'x');
% set(groot,'defaultLineLineWidth',defaultLineWidth);

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