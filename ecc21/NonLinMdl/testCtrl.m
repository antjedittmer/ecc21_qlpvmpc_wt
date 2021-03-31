%% Clear workspace
clear; close all; clc; 
restoredefaultpath;
bdclose all;
clear ABgL_WT; % LUT variables saved as persistent variables
Simulink.data.dictionary.closeAll;

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
subDir = {'subfunctionsMPC', 'subfunctionsUtil', 'WECSCtrl', 'WECSMDl','WECSAct'}; 

% Get relative paths with respect to testCrl.m file
fileName = mfilename('fullpath');
workDir = fileparts(fileName);
parentDir = fileparts(workDir);

% Add parentDir and workdir to path: Can be called in both directories
addpath(parentDir)
addpath(workDir)

% Add LinMdl (for WECS param in workspace for testing)
addpath(fullfile(parentDir,'LinMdl'));

% Add path to input data (extracted from FAST with FASTtool)
dataInPath = fullfile(parentDir,dataIn);
addpath(dataInPath);

% Add output path to save mat files
dataOutPath = fullfile(workDir,dataOut); 
if ~isdir(dataOutPath) %#ok<ISDIR>
    mkdir(dataOutPath);
end
addpath(dataOutPath);

% Add subdirerectoreis
for idx = 1: length(subDir)
    addpath(fullfile(workDir,subDir{idx}));
end

%% Load FASTTool wind as test signal
% onlyInit = 1; DT = 0.008; %Nout used anymore
%load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');
load(fullfile(dataInPath,'OutDataWind18NTW.mat'),'OutTable');


