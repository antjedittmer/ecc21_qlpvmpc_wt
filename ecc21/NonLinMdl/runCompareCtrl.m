%clear ; close all; clc;
function runCompareCtrl(loadData,strWindType)



% Set path to inputdata directory
%% Handle optional inputs
if ~nargin || isempty(loadData) 
   loadData = 1; %load simulation output data if available;
end

% Two testcases: step sweep 4 to 25 ms and normal dist. 18 m/s mean
if nargin < 2 || isempty(strWindType)
    strWindType = 'Sweep'; %'NTW18';
end

% Set path to data output directory
%loadData = 1; %load simulation output data if available;
dataDirOut = 'dataOut';
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Load data closed loop simulation PI


% Weight values (copied from MATLAB function
q_ = [1 10^4 0 10^3 10^3 0 0 0]; 
r_ = [1 10^4]; %[10^7,0.1];
p = 10^3;

%% 
% Provide names of Simulink models to be run
if strcmp(strWindType,'Sweep') == 1
    % sweep from 4 to 25 in steps
    outDataSimulationMat = 'OutDataSweep.mat';
    strFig = '';
    testCaseStr = 'Wind Sweep';
else
    % wind with average 18 m/s
    outDataSimulationMat = 'OutDataWind18NTW.mat';
    strFig = 'NTW18';
    testCaseStr = 'Wind, mean 18 m/s';
end


load(outDataSimulationMat ,'OutTable');%

OutTableTest2 = OutTable; %(:,OutTableTest2.Properties.VariableNames(2:end));
OutList = OutTableTest2.Properties.VariableNames;
idxWind = contains(OutList, 'Wind');
OutTableTest2.Wind = sqrt(sum((OutTable{:,idxWind}.^2),2));
OutTableTest2.GenTq = OutTableTest2.GenTq*1000; % kNm -> Nm
OutTableTest2.RotSpeed = OutTableTest2.RotSpeed *pi/30; %RPM -> rad/s

%%
simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta.slx'; %
if strcmp(strWindType,'Sweep') == 1
    
    % Names of simulation output in correct order
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', ...
    'NcIMUTAxs', 'NcIMUTAys'}; %%'Tower FA acc', 'Tower SS acc'};

% Run simulation for two models or load mat files if available

matFileOutTableTest1 = fullfile(dataDirOut,['OutTableMPC',strFig,'.mat']);
OutTableTest1 = getSimulationOutputTable(matFileOutTableTest1,loadData,OutTableTest2,simMdlname,varnames);

    load('OutDataTestCL_NoFilt_Simulink_beta_steadyStateCLSim.mat',...
        'OutDataTestOL','tictoc_LPVMPC');
    
    figStr = 'RefSimulinkClosedLoop_beta';
    meanOpt = strrep(sprintf(' mean cputime MPC: %1.2e s',mean(tictoc_LPVMPC.Data)),'e-0','e-');
    k_pitch = 1;
    timeForPlot = [15,400];
    plotComparison(OutTableTest2,OutDataTestOL,q_,r_,p,DT,figStr,figDir,rOut,timeForPlot,k_pitch,'',meanOpt)
    
else
    load('OutDataTestCLNTW18_NoFilt_Simulink_beta_steadyStateCLSim.mat','OutDataTestOL','q_','r_','p','rOut', 'tictoc_LPVMPC');
    figStr = 'RefSimulinkClosedLoopNTW18_beta';
    k_pitch = 1;
    timeForPlot = [15,400];
    meanOpt = strrep(sprintf(' mean cputime MPC: %1.2e s',mean(tictoc_LPVMPC.Data)),'e-0','e-');
    plotComparison(OutTableTest2,OutDataTestOL,q_,r_,p,DT,figStr,figDir,rOut,timeForPlot,k_pitch,'',meanOpt);
    
end


function plotComparison(OutTableTest2,OutDataTestOL,q_,r_,p,DT,figStr,figDir,r,timeForPlot,k_pitch,plot6,meanOpt)

if nargin < 9
    r = [];
end

maxTime = min(height(OutTableTest2),length(OutDataTestOL));
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

varnames = OutTableTest2.Properties.VariableNames;
OutTableTestCL = array2table(OutDataTestOL,'VariableNames',varnames);

idxPlot = timeVec >= timeForPlot(1) & timeVec <= timeForPlot(2);
timeVecPlot = timeVec(idxPlot);
OutTableTest2Plot = OutTableTest2(idxPlot,:);
OutTableTestCLPlot = OutTableTestCL(idxPlot,:);

% Calculate omega_g
Generator_Efficiency = 0.944;
omegaG1 = 1000 /Generator_Efficiency * OutTableTest2Plot.GenPwr./OutTableTest2Plot.GenTq;
omegaG2 = 1000 /Generator_Efficiency * OutTableTestCLPlot.GenPwr./OutTableTestCLPlot.GenTq;

% Set reference
idxVec = 20/DT : 20/DT : length(timeVecPlot)-50;
LUT_wind2 = round(OutTableTest2Plot.Wind(idxVec));
LUT_torque2 = OutTableTest2Plot.GenTq(idxVec); % kNm (in Simulink
LUT_omegaG2 = OutTableTest2Plot.RotSpeed(idxVec) * pi/30 *97; %RPM -> rad/s
%LUT_pitch2 = OutTableTest2Plot.BlPitch1(idxVec)* pi/180;

if ~isempty(r)
    TGRef = r(idxPlot,1);
    try
        omegaGRef = interp1(LUT_wind2,LUT_omegaG2, OutTableTest2Plot.Wind,'spline');
        
    catch
        omegaGRef = mean(omegaG2)*ones(size(OutTableTest2Plot.Wind));% LUT_omegaG2(end)* ones(size(OutTableTest2Plot.Wind));
        
    end
    
    
else
    try
        omegaGRef = interp1(LUT_wind2,LUT_omegaG2, OutTableTest2Plot.Wind,'spline');
        TGRef = interp1(LUT_wind2,LUT_torque2, OutTableTest2Plot.Wind,'spline');
    catch
        omegaGRef = LUT_omegaG2(end)* ones(size(OutTableTest2Plot.Wind));
        TGRef = LUT_torque2(end)* ones(size(OutTableTest2Plot.Wind));
    end
end
PGRef = TGRef.*omegaGRef* 0.944;

% Strings of weights for optimization
qVec = deblank(sprintf('%d ',ceil(q_)));
rVec = deblank(sprintf('%d ', ceil(r_)));
pVec = deblank(sprintf('%d', p(1)/q_(1)));

% yAxCell = {'V [m/s]', '\omega_g [rad/s]', 'T_g [kNm]', 'Pitch \beta [°]',...
%     'P_g [kW]','y_t [m/s^2]'};

yAxCell = {'wind V [m/s]', 'RotSpd \omega_g [rad/s]', 'GenTq T_g [kNm]', 'Pitch \beta [°]',...
    'GenPwr P_g [kW]','Twr_{FA} y_t [m/s^2]'};

figure;
if isempty(plot6)
    ax1(1) = subplot(3,1,1);
else
    ax1(1) = subplot(6,1,1);
end
plot(timeVecPlot, OutTableTest2Plot.Wind);
ylabel(yAxCell{1}); %'V [m/s]')
axis tight; grid on;

varPI = var(OutTableTest2Plot.NcIMUTAxs);
varqLPV = var(OutTableTestCLPlot.NcIMUTAxs);
titleStdPI = sprintf('%2.2e',varPI);
titleStdqLPV = sprintf('%2.2e',varqLPV);
ratioStd = sprintf('%2.1f ', varPI/varqLPV);

varPI_Pwr = var(abs(PGRef/1000 - OutTableTest2Plot.GenPwr));
varqLPV_Pwr = var(abs(PGRef/1000 - OutTableTestCLPlot.GenPwr));
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
    plot(timeVecPlot, OutTableTest2Plot.GenTq/1000, timeVecPlot,OutTableTestCLPlot.GenTq/1000,'-.',...
        timeVecPlot, TGRef/1000, 'k--');
    posAxis = axis;
    axis([posAxis(1:2), min(40,posAxis(3)), 45]);
    ylabel(yAxCell{3}); %'T_g [kNm]')
    grid on; % axis tight;
    
    ax1(4) = subplot(6,1,4);
    plot(timeVecPlot, OutTableTest2Plot.BlPitch1, timeVecPlot,OutTableTestCLPlot.BlPitch1 * k_pitch,'-.');
    ylabel(yAxCell{4}); %'Pitch \beta [°]')
    axis tight; grid on;
end

if isempty(plot6)
    ax1(2) = subplot(3,1,2);
else
    ax1(5) = subplot(6,1,5);
end
plot(timeVecPlot, OutTableTest2Plot.NcIMUTAxs, timeVecPlot,OutTableTestCLPlot.NcIMUTAxs,'-.');
title([yAxCell{6}, ': var_{P-Pi} = ',titleStdPI , ', var_{qLPVMPC} = ' ,titleStdqLPV,', ratio: ' ,ratioStd])
ylabel(yAxCell{6}); %'y_t [m/s^2]')
axis tight; grid on;


if isempty(plot6)
    ax1(3) = subplot(3,1,3);
else
    ax1(6) = subplot(6,1,6);
end
plot(timeVecPlot, OutTableTest2Plot.GenPwr, timeVecPlot,OutTableTestCLPlot.GenPwr,'-.',timeVecPlot,PGRef/1000,'k:');
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

%snapnow;
print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),plot6,'_',figStr,'_',num2str(timeVecPlot(end))], '-dpng');
print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),plot6,'_',figStr,'_',num2str(timeVecPlot(end))], '-depsc');

