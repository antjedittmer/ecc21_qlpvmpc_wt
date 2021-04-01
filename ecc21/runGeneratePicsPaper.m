%% Script to generate pictures for paper 

%% Set path to initialization script testCtrl and run it
addpath('NonLinMdl');
testCtrl;

%% Matlab analysis: Bode plots linearized FAST models (reference)
compareLinearModels;

%% Simulink simulations 

% Run models in closed loop w baseline controller (k-omega-squared - Pi)
close all;
runCompareModels('NTW18',1);
close all;
runCompareModels('Sweep',1);

% Run models in closed loop w qLPV MPC
close all;
runCompareCtrl('NTW18',1);

close all;
runCompareCtrl('Sweep',1);
