%% Clear workspace
clear; close all; clc; 
restoredefaultpath;
bdclose all;
% Simulink.data.dictionary.closeAll;

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
dirSub = {'subfunctionsUtil', 'subfunctionsMPC','WECSCtrl', 'WECSMDl','WECSAct'}; 

% Get relative paths with respect to testCrl.m file
dirFileName = mfilename('fullpath');
dirWork = fileparts(dirFileName);
dirParent = fileparts(dirWork);

% Add parentDir and workdir to path: Can be called in both directories
addpath(dirParent)
addpath(dirWork)

% Add LinMdl (for WECS param in workspace for testing)
addpath(fullfile(dirParent,'LinMdl'));

% Add path to input data (extracted from FAST with FASTtool)
dataInPath = fullfile(dirParent,dataIn);
addpath(dataInPath);

% Add output path to save mat files
dataOutPath = fullfile(dirWork,dataOut); 
if ~isdir(dataOutPath) %#ok<ISDIR>
    mkdir(dataOutPath);
end
addpath(dataOutPath);

% Add subdirerectoreis
for dirIdx = 1: length(dirSub)
    addpath(fullfile(dirWork,dirSub{dirIdx}));
end

%% Load FASTTool wind as test signal
DT = 0.008; % one sample time 
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');
%load(fullfile(dataInPath,'OutDataWind18NTW.mat'),'OutTable');

%% Create bus objects used in all Simulink models
% These buses do not need to be created anymore as they are stored in the
% data dictionaries
Mdl_TuDelftOLFunction(0);
MdLCtrl_BaselineFunction(0);

%% Initialize path to FASTFile: FAST_InputFileName and OutList
testqLPV; %
TMax = 660; %Can be overwritten later on


%% Clean-up
clear RESTORE*
%clear data*
clear dir*



