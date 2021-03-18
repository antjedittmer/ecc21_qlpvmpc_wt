%% get values
clear;
clear ABgL_WT;
close all;
Simulink.data.dictionary.closeAll

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
subDir = {'subfunctionsMPC', 'WECSCtrl', 'WECSMDl'}; 

dataInPath = fullfile(fileparts(pwd),dataIn);
addpath(dataInPath);

if ~isdir(dataOut) %#ok<ISDIR>
    mkdir(dataOut);
end
addpath(fullfile(pwd,dataOut));

for idx = 1: length(subDir)
    addpath(subDir{idx});
end

onlyInit = 1;
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable'); %#ok<NASGU>
%DT = 0.008;

