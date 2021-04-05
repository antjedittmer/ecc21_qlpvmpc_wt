function runCompareCtrl(strWindType,loadData,figNo1,useFASTforComparison)
% runComparCtrl compares simulation results for the baseline and qLPV MPC
% controller in closed loop with the Simulink WECS model.
%
% All inputs are optional:
% - strWindType: Two testcases: step sweep 4 to 25 ms and normal dist. with 
%   18 m/s mean (Default: Sweep)
% - loadData: load simulation output data if available instead of running
%   simulation (Default: 1)
% - figNo1: Number of figure (Default: 1)
% - useFASTforComparison: Uses FAST(1) or Simulink data for comparison

%% Handle optional inputs
% The default inputs are provided here.

if ~nargin || isempty(strWindType) 
    strWindType = 'Sweep'; %'NTW18';   
end

% Two testcases: step sweep 4 to 25 ms and normal dist. 18 m/s mean
if nargin < 2 || isempty(loadData)
    loadData = 1; %load simulation output data if available;
end

if nargin < 3 || isempty(figNo1) 
    figNo1 = 1;
end

if nargin < 4 || isempty(useFASTforComparison)
    useFASTforComparison = 0;
end

%% Set path to figure directory

% Set path to inputdata directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);
dataInDir = fullfile( mainDir,'dataIn');
addpath(dataInDir);

% Set path to figure directory
figDir = fullfile(mainDir,'figDir');
if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory
dataDirOut = fullfile(workDir, 'dataOut');
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Names of simulation output in correct order
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', ...
    'NcIMUTAxs', 'NcIMUTAys'}; %%'Tower FA acc', 'Tower SS acc'};

% Weight values (copied from MATLAB function)
q_ = [1 10^4 0 10^3 10^3 0 0 0];
r_ = [1 10^4]; %[10^7,0.1];
p = 10^3;

% Switch between Sweep and NTW18
if strcmp(strWindType,'Sweep') == 1 % sweep from 4 to 25 in steps
    outDataSimulationMat = 'OutDataSweep.mat';
    strFig = ''; %testCaseStr = 'Wind Sweep';
else % wind with average 18 m/s
    outDataSimulationMat = 'OutDataWind18NTW.mat';
    strFig = 'NTW18'; %testCaseStr = 'Wind, mean 18 m/s';
end

%% Load baseline and MPC data or run simulations

%For FAST simulation
load(outDataSimulationMat ,'OutTable'); 

% Load data closed loop simulation PI
if useFASTforComparison % from FAST
    OutTableTest2 = OutTable; 
else % from Simulink
    loadBaseline = 1; %always load if available
    matFileOutTableTest2 = fullfile(dataDirOut,['OutTableTest2',strFig,'.mat']);
    simMdlname2 = 'test_SimulinkMdl2_Baseline';
    OutTableTest2 = getSimulationOutputTable(matFileOutTableTest2,loadBaseline,OutTable,simMdlname2);
end

% Check that SI units are used (can be removed, 
if mean(OutTableTest2.GenTq) <100 % protection against legacy data in kNm
    OutTableTest2.GenTq = OutTableTest2.GenTq*1000; % kNm -> Nm
end

if mean(OutTableTest2.RotSpeed) >1 
    OutTableTest2.RotSpeed = OutTableTest2.RotSpeed *pi/30; %RPM -> rad/s
end

% Run or load qLPV controller
simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta.slx'; %
matFileOutTableTest1 = fullfile(dataDirOut,['OutTableMPC',strFig,'.mat']);
[OutTableMPC, tictoc_LPVMPC,GenPwrRef] = getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname,varnames);

%% Prepare information for plot and plot data

% Information for plot
figStr = ['RefSimulinkClosedLoop',strFig,'_beta'];
meanOpt = strrep(sprintf(' mean cputime MPC: %1.2e s',mean(tictoc_LPVMPC.Data)),'e-0','e-');
k_pitch = 1;
timeForPlot = [15,400];
DT = 0.008;
r = GenPwrRef;

% Plot
plotComparison(OutTableTest2,OutTableMPC,q_,r_,p,DT,figStr,figDir,r,timeForPlot,k_pitch,'', meanOpt,figNo1)


function plotComparison(OutTableTest2,OutTableMPC,q_,r_,p,DT,figStr,figDir,r,timeForPlot,k_pitch,plot6,meanOpt,figNo1)

if nargin < 9
    r = [];
end

maxTime = min(height(OutTableTest2),height(OutTableMPC));
if nargin <10
    timeForPlot(1) = 0;
    timeForPlot(2) = maxTime;
end

if nargin <11
    k_pitch = 180/pi;
end

if nargin <12
    plot6 = '';
end

if nargin < 13
    meanOpt = '';
end

timeVec = 0:DT: maxTime*DT - DT;

idxPlot = timeVec >= timeForPlot(1) & timeVec <= timeForPlot(2);
timeVecPlot = timeVec(idxPlot);
OutTableTest2Plot = OutTableTest2(idxPlot,:);
OutTableMPCPlot = OutTableMPC(idxPlot,:);

% Calculate omega_g
Generator_Efficiency = 0.944;
omegaG1 = 1000 /Generator_Efficiency * OutTableTest2Plot.GenPwr./OutTableTest2Plot.GenTq;
omegaG2 = 1000 /Generator_Efficiency * OutTableMPCPlot.GenPwr./OutTableMPCPlot.GenTq;

PGRef = r.Data(idxPlot);

% Strings of weights for optimization
qVec = deblank(sprintf('%d ',ceil(q_)));
rVec = deblank(sprintf('%d ', ceil(r_)));
pVec = deblank(sprintf('%d', p(1)/q_(1)));

% yAxCell = {'V [m/s]', '\omega_g [rad/s]', 'T_g [kNm]', 'Pitch \beta [°]', 'P_g [kW]','y_t [m/s^2]'};
yAxCell = {'wind V [m/s]', 'RotSpd \omega_g [rad/s]', 'GenTq T_g [kNm]', 'Pitch \beta [°]',...
    'GenPwr P_g [kW]','Twr_{FA} y_t [m/s^2]'};

figure(figNo1);
if isempty(plot6)
    ax1(1) = subplot(3,1,1);
else
    ax1(1) = subplot(6,1,1);
end
plot(timeVecPlot, OutTableMPCPlot.Wind);
ylabel(yAxCell{1}); %'V [m/s]')
axis tight; grid on;

varPI = var(OutTableTest2Plot.NcIMUTAxs);
varqLPV = var(OutTableMPCPlot.NcIMUTAxs);
titleStdPI = sprintf('%2.2e',varPI);
titleStdqLPV = sprintf('%2.2e',varqLPV);
ratioStd = sprintf('%2.1f ', varPI/varqLPV);

varPI_Pwr = var(abs(PGRef - OutTableTest2Plot.GenPwr));
varqLPV_Pwr = var(abs(PGRef - OutTableMPCPlot.GenPwr));
titleStdPI_Pwr = sprintf('%2.2e',varPI_Pwr);
titleStdqLPV_Pwr = sprintf('%2.2e',varqLPV_Pwr);
ratioStd_Pwr = sprintf('%2.1f ', varPI_Pwr/varqLPV_Pwr);

title({['Baseline (P-PI) vs. qLPV MPC; ', meanOpt],... %WT mdl: Model2
    ['q_ = [', qVec,'], r = [',rVec,'], p = ',pVec,'* q']})

if ~isempty(plot6)
    ax1(2) = subplot(6,1,2);
    plot(timeVecPlot, omegaG1,timeVecPlot, omegaG2,'-.', timeVecPlot, omegaGRef,'k--');
    ylabel(yAxCell{2}); %'\omega_g [rad/s]')
    axis tight; grid on;
    
    ax1(3) = subplot(6,1,3);
    plot(timeVecPlot, OutTableTest2Plot.GenTq/1000, timeVecPlot,OutTableMPCPlot.GenTq/1000,'-.',...
        timeVecPlot, TGRef/1000, 'k--');
    posAxis = axis;
    axis([posAxis(1:2), min(40,posAxis(3)), 45]);
    ylabel(yAxCell{3}); %'T_g [kNm]')
    grid on; % axis tight;
    
    ax1(4) = subplot(6,1,4);
    plot(timeVecPlot, OutTableTest2Plot.BlPitch1, timeVecPlot,OutTableMPCPlot.BlPitch1 * k_pitch,'-.');
    ylabel(yAxCell{4}); %'Pitch \beta [°]')
    axis tight; grid on;
end

if isempty(plot6)
    ax1(2) = subplot(3,1,2);
else
    ax1(5) = subplot(6,1,5);
end
plot(timeVecPlot, OutTableTest2Plot.NcIMUTAxs, timeVecPlot,OutTableMPCPlot.NcIMUTAxs,'-.');
title([yAxCell{6}, ': var_{P-Pi} = ',titleStdPI , ', var_{qLPVMPC} = ' ,titleStdqLPV,', ratio: ' ,ratioStd])
ylabel(yAxCell{6}); %'y_t [m/s^2]')
axis tight; grid on;

if isempty(plot6)
    ax1(3) = subplot(3,1,3);
else
    ax1(6) = subplot(6,1,6);
end
plot(timeVecPlot, OutTableTest2Plot.GenPwr, timeVecPlot,OutTableMPCPlot.GenPwr,'-.',timeVecPlot,PGRef,'k:');
title(['|P_{g,ref}  - P_g| [kW]: var_{P-Pi} = ',titleStdPI_Pwr, ', var_{qLPVMPC} = ' ,...
    titleStdqLPV_Pwr,', ratio: ' ,ratioStd_Pwr])
ylabel(yAxCell{5}); %'P_g [kW]')
axis tight; grid on;
legend('P-PI','MPC','P_{g,ref}','Location','SouthEast')

xlabel('time [s]');

linkaxes(ax1,'x')
set(gcf,'Name', figStr)

posDefault = get(gcf, 'position');
if ~isempty(plot6)
    set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.5]);
else
    set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.2]);
end

print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),plot6,'_',figStr,'_',num2str(timeVecPlot(end))], '-dpng');
print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),plot6,'_',figStr,'_',num2str(timeVecPlot(end))], '-depsc');

