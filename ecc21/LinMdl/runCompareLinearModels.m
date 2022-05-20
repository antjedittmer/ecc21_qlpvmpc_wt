
close all; clear; clc;
speedVec = [4,5,10,11,12,15,22]; % %18-3;

useActuatorStates = 0;
figNoAdd = 0;
createBodePlots = 1;
plotVisible = 1;
noOut = 3;

figFolder = 'figLinUIO1';

%compareLinearModels(speedVec,figFolder, useActuatorStates,figNoAdd,createBodePlots,plotVisible)
sysOut = compareLinearModels(speedVec,figFolder,useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut);