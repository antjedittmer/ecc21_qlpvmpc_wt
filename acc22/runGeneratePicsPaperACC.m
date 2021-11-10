%% Script to generate pictures for paper
clc; clear; close all;
restoredefaultpath;

%% Set path to initialization script initWorkspace and run it
addpath('LinMdl');

%% Create cq and ct curves
plotOn = 1;
initModel5MWNREL(plotOn);

%% 
addpath('NonLinMdl');
%initWorkspace;
initUIO_WECS;


%% Simulink simulations 
% Simulink models are compared with FAST (NREL) references.
% FAST simulation results obtained with FASTTool (Tu Delft) are provided as
% mat-files in dataIn folder. 

%Load data if available from previous simulation.
strWindType = ''; %Default wind: nomal turbulent wind, 18 
loadData = 1; %use data if available
useFASTForComparison = 1; % use data recorded with FAST in closed loop with PI

figNo1 = 100;
UIOstr = '';
runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,UIOstr);


figNo1 = 101;
UIOstr = 'UIO';
runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,UIOstr);




