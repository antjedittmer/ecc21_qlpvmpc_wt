
mainDir = fullfile(pwd);
addpath(fullfile(fileparts(pwd),'dataIn'));
addpath(fullfile(pwd,'dataOut'));

% Set path to figure directory
%figDir = fullfile(mainDir,'html','figDir');
figDir = fullfile(fileparts(mainDir),'figDir');
if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory
%loadData = 1; %load simulation output data if available;
dataDirOut = 'dataOut';
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

addpath('subfunctions');

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
load(matFileOutTableTest1,'OutTableTest1');

load('OutDataTestCLNTW18_NoFilt_Simulink_beta_steadyStateCLSim.mat','OutDataTestOL','q_','r_','p');
figStr = 'RefSimulinkClosedLoopNTW18_beta';
k_pitch = 1;

maxTime = min(height(OutTableTest1),length(OutDataTestOL));

DT = 0.008;
timeVec = 0:DT: maxTime*DT - DT;
OutTableTest1.Time = timeVec';
OutTableTestCL.Time = timeVec';

cellOfTable{1} = OutTableTest1;
cellOfTable{2} = OutTableTest1; %OutTableTestCL;
calculateEvalCriteria(cellOfTable)


