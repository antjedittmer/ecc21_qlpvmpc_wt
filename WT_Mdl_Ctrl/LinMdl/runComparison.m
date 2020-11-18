speedVec = [1,10,22];
useTUDelftParam = 1;
%figFolder = 'FigureLagrange_201015_comp';
useActuatorStates = 0;

figFolder = fullfile(fileparts(mainDir),'figDir');

figNoAdd = 0;
Linear_Models_Lagrange_Comparison(speedVec,useTUDelftParam,figFolder, useActuatorStates,figNoAdd)

% figNoAdd = 1;
% useTUDelftParam = 0;
% Linear_Models_Lagrange_Comparison(speedVec,useTUDelftParam,figFolder, useActuatorStates,figNoAdd)