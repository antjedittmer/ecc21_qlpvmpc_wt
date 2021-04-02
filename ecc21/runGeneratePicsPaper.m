%% Script to generate pictures for paper 

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
% type of wind (either a step sweep from 4 to 25 ms or normal turbulent
% wind with mean 18 ms)

% Run Simulink models in closed loop w baseline controller( Torque controller
% k-omega-squared, Pitch controller: Gainscheduled Pi)
close all;
runCompareModels('NTW18',1);
close all;
runCompareModels('Sweep',1);

% Run models in closed loop w qLPV MPC
close all;
runCompareCtrl('NTW18',1);

close all;
runCompareCtrl('Sweep',1);
