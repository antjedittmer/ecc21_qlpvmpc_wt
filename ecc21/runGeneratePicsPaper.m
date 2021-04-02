%% Script to generate pictures for paper
clc; clear; close all;

%% Set path to initialization script testCtrl and run it
addpath('NonLinMdl');
testCtrl;

%% Matlab analysis: Bode plots linearized FAST models (reference)
compareLinearModels;

%% Simulink simulations 
% Simulations are run to compare the Simulink models with FAST reference.
% FAST simulation results obtained with FASTTool (Tu Delft) are provided as
% mat-files in dataIn folder. 
% Inputs to the files runCompareModels and runCompareCtrl are 
% - strWindType: Two testcases: step sweep 4 to 25 ms and normal dist. with 
%   18 m/s mean (Default: Sweep)
% - loadData: load simulation output data if available instead of running
%   simulation (Default: 1)
loadData = 1;
% Run Simulink models in closed loop w baseline controller( Torque controller
% k-omega-squared, Pitch controller: Gainscheduled Pi)
figNo = 1;
runCompareModels('Sweep',loadData,figNo);
figNo = 4;
runCompareModels('NTW18',loadData,figNo);


% Run models in closed loop w qLPV MPC
figNo1 = 7;
runCompareCtrl('Sweep',loadData,figNo1);
figNo1 = 8;
runCompareCtrl('NTW18',loadData,figNo1);


