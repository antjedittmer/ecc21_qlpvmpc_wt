
%% Add all relevant paths

% Matlab functions
addpath('subfunctions');

mainDir = fullfile(pwd);
addpath(fullfile(fileparts(pwd),'dataIn'));
dataDirOut = 'dataOut1';


% Set path to figure directory
%figDir = fullfile(mainDir,'html','figDir');
figDir = fullfile(fileparts(mainDir),'figDir');
if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory
%loadData = 1; %load simulation output data if available;
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end
addpath(fullfile(pwd,dataDirOut));

% Load data closed loop simulation PI
DT = 0.008;
strWindType = 'NTW';
if strcmp(strWindType,'Sweep') == 1
    % sweep from 4 to 25 in steps
    workspaceSimulationMat = 'workSpace.mat';
    outDataSimulationMat = 'OutDataSweep.mat';
    strFig = '';
    testCaseStr = 'Wind Sweep';
else
    % wind with average 18 m/s
    workspaceSimulationMat = 'workSpace18NTM.mat';
    outDataSimulationMat = 'OutDataWind18NTM.mat';
    strFig = 'NTW18';
    testCaseStr = 'Wind, mean 18 m/s';
end
matFileOutTableTest1 = fullfile(dataDirOut,['OutTableTest1',strFig,'.mat']);
load(matFileOutTableTest1,'OutTableTest1NTW18');
OutTableTest1=  OutTableTest1NTW18;

load('OutDataTestCLNTW18_NoFilt_Simulink_beta_steadyStateCLSim.mat','OutDataTestOL','q_','r_','p');
figStr = 'RefSimulinkClosedLoopNTW18_beta';
k_pitch = 1;

maxTime = min(height(OutTableTest1),length(OutDataTestOL));

timeVec = 0:DT: maxTime*DT - DT;
OutTableTest1.Time = timeVec';
OutTableTestCL.Time = timeVec';

cellOfTable{1} = OutTable; %Test1;
cellOfTable{2} = OutTable; %Test1; %OutTableTestCL;

tmp = load('cellOfTable');
x = load('OutTableTest_rate13_IPC.mat');
load('OutTableTest_rate08_PI.mat','OutTableTest_rate08_PI');
load('OutTableTest_rate01_MPC.mat','OutTableTest_rate01_MPC');
load('OutTableTest_rate13_MPC.mat','OutTableTest_rate13_MPC');

cellOfTable4{1} = cellOfTable{1};
cellOfTable4{2} = cellOfTable{3};
strCell = {'PI-CPC','PI-IPC'}; %,'MPC, 13deg/s','MPC, 1deg/s'};
calculateEvalCriteria(cellOfTable4,strCell)

OutTable4 = OutTableTest_rate08_PI;
OutTable4.Time = cellOfTable{1}.Time(1:height(OutTable4));

cellOfTable5{1} = OutTableTest_rate08_PI{1};
cellOfTable5{2} = cellOfTable{3};


