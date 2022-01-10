%% Script to generate pictures for paper
clc; clear; close all;
restoredefaultpath;

%% Set path to initialization script initWorkspace and run it
addpath('LinMdl');

%% Create cq and ct curves
plotOn = 1;
initModel5MWNREL(plotOn);

%% 
addpath(genpath('NonLinMdl'));
%initWorkspace;
initUIO_WECS;


%% Simulink simulations 
% Simulink models are compared with FAST (NREL) references.
% FAST simulation results obtained with FASTTool (Tu Delft) are provided as
% mat-files in dataIn folder. 

%Load data if available from previous simulation.
strWindType = ''; %Default wind: nomal turbulent wind, 18 
useFASTForComparison = 1; % use data recorded with FAST in closed loop with PI

figNo1 = 100;
UIOstr = '';
loadData = 1;
runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,UIOstr);

figNo1 = 101;
switchNoise = 0;% run Noise with limit 10^-6
UIOstr = 'UIO';
runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,UIOstr);

figNo1 = 103;
UIOstr = 'UIOall';
[OutTableTest2Plot, OutTableMPCPlot,OutTableMPC2Plot,PGRef] = runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,UIOstr);



figure(104) 
%workDir = fileparts(mfilename('fullpath'));
mainDir = pwd; %fileparts(workDir);

% Set path to figure directory
figDir = fullfile(mainDir,'figDirUIO');
subplot(1,2,1);
histogram(OutTableTest2Plot.NcIMUTAxs,20); hold on;
histogram(OutTableMPCPlot.NcIMUTAxs,20); histogram(OutTableMPC2Plot.NcIMUTAxs,20); 
hold off;
axis tight; grid on;
ylabel('Occurence [-]');
xlabel( 'Twr_{FA} y_t [m/s^2]')
subplot(1,2,2);
histogram(abs(OutTableTest2Plot.GenPwr - PGRef),20); hold on;
histogram(abs(OutTableMPCPlot.GenPwr - PGRef),20); histogram(abs(OutTableMPC2Plot.GenPwr - PGRef),20); 
xlabel('|P_{g,ref}-P_g| [kW]')
% histogram(OutTableTest2Plot.GenPwr - PGRef,20); hold on;
% histogram(OutTableMPCPlot.GenPwr - PGRef,20); histogram(OutTableMPC2Plot.GenPwr - PGRef,20); 
% xlabel('P_{g,ref}-P_g [kW]')
legend('P-PI','MPC V_{est}','MPC V_{real}','Position',[0.3285 0.7738 0.1792 0.1341])
hold off;
axis tight; grid on;
% save to png/eps
print(gcf,fullfile(figDir,'histTwrPabs'), '-dpng');
print(gcf,fullfile(figDir,'histTwrPabs'), '-depsc');






