addpath(fullfile(fileparts(pwd),'dataIn'));
initModel5MWNRELBianchi;
load('workSpace.mat')
load('OutDataSweep.mat');
OutTable = array2table(OutData,'VariableNames', OutList);
idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));