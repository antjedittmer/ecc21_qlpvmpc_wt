%% get values
clear;
clear ABgL_WT;
close all;
Simulink.data.dictionary.closeAll;

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
subDir = {'subfunctionsMPC', 'WECSCtrl', 'WECSMDl'}; 


fileName = mfilename('fullpath');
workDir = fileparts(fileName);
parentDir = fileparts(workDir);

dataInPath = fullfile(parentDir,dataIn);
addpath(dataInPath);

dataOutPath = fullfile(workDir,dataOut); 
if ~isdir(dataOutPath) %#ok<ISDIR>
    mkdir(dataOutPath);
end
addpath(dataOutPath);

for idx = 1: length(subDir)
    addpath(fullfile(workDir,subDir{idx}));
end

onlyInit = 1;
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable'); %#ok<NASGU>
%DT = 0.008;

