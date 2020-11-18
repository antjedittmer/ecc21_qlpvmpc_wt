speedVec = [1,10,22];
useTUDelftParam = 1;
%figFolder = 'FigureLagrange_201015_comp';
useActuatorStates = 0;
mainDir = fullfile(pwd);
dataInDir = fullfile(fileparts(pwd),'dataIn');

figDir = fullfile(fileparts(mainDir),'figDir');
if ~isfolder(figDir)
    mkdir(figDir)
end

figNoAdd = 0;
Linear_Models_Lagrange_Comparison(speedVec,useTUDelftParam,figDir, useActuatorStates,figNoAdd)

% figNoAdd = 1;
% useTUDelftParam = 0;
% Linear_Models_Lagrange_Comparison(speedVec,useTUDelftParam,figFolder, useActuatorStates,figNoAdd)