function plotOutTable(OutTable,OutTableTest1,OutTableTest2,testCaseStr, testCaseCell)

yAxCell = {'wind V [m/s]', 'GenTq T_g [kNm]', 'Pitch \beta [°]', 'RotSpd \omega_r [rpm]',...
    'GenPwr P_g [MW]','Twr_{FA} y_t [m/s^2]'};
if nargin < 5 || isempty(testCaseCell)
    testCaseCell = {'Mdl1: Rot+Twr ', 'Mdl2: Rot,Twr,Bld+Act ','FAST'}
end

figNo2 = 1;
strFig = 'testConstr';
% Assign data from FAST to variables
idxTime = 1: height(OutTable);
OutTable = OutTable(idxTime,:);
OutTable.BldPitch1 = OutTable.BlPitch1;
OutTable.Torque = OutTable.GenTq;

% Calculate and assign wind amplitude from comonents
idxWind = contains(OutTable.Properties.VariableNames, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

% Get color map for 'title legends'
cl = colormap('lines');
% titleStr2 = [testCaseStr,': ',testCaseCell{3},' {\color[rgb]{',num2str(cl(1,:)),'}',testCaseCell{3},' ',...
%     '\color[rgb]{',num2str(cl(2,:)),'}',testCaseCell{1}];
% titleStr = [testCaseStr,': Mdl1: Rot+Twr {\color[rgb]{',num2str(cl(1,:)),'}Mdl2: Rot,Twr,Bld+Act ',...
%     '\color[rgb]{',num2str(cl(2,:)),'}FAST} '];
titleStr = [testCaseStr,': ',testCaseCell{1},'{\color[rgb]{',num2str(cl(1,:)),'}',testCaseCell{2},...
    '\color[rgb]{',num2str(cl(2,:)),'}',testCaseCell{3},'}'];

% Create output plot time reference
idxTime = 1: min([length(idxTime), height(OutTableTest1),height(OutTableTest2)]);
time = OutTable.Time(idxTime);
idxTime = time > 15;
time = time(idxTime);


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

set(gcf,'Name',['cmpTimeDomain_All',strFig])
posDefault = get(gcf, 'position');
% set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.5]);
% set(groot,'defaultLineLineWidth',defaultLineWidth);
% 
% print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-dpng');
% print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');