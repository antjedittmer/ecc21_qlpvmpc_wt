%% Clear workspace
clear; close all; clc; 
clear ABgL_WT; % LUT variables saved as persistent variables
Simulink.data.dictionary.closeAll;

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
subDir = {'subfunctionsMPC', 'WECSCtrl', 'WECSMDl','WECSAct'}; 

% Get relative paths with respect to testCrl.m file
fileName = mfilename('fullpath');
workDir = fileparts(fileName);
parentDir = fileparts(workDir);

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
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');


