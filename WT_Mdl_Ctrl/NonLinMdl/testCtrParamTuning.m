clc; clear; bdclose;

% initialize constants to base workspace
testCtrl;
load(fullfile(dataInPath,'OutDataWind18NTM.mat'),'OutData','OutTable');

onIPC = 0;
switchIPC = 1;
initDLRparameter(onIPC,switchIPC);

[sysd,x0,uoff,xoff,yoff] = initLPVParam;

sim('test_SimulinkMdl3LPV_IPCa');
OutTable = array2table(OutDataLPV,'VariableNames',OutList);

%idxRootMyb = contains(OutList,'RootMyb') | contains(OutList,'RootMFlp');
%OutTable{:,idxRootMyb} = MY;

[W,K,ref] = getEvalCriteria(OutTable);
% Put(strmatch('RootMFlp1',OutList))OutDataLPV()

onIPC = 1;

switchIPC = 1;
initDLRparameter(onIPC,switchIPC);
sim('test_SimulinkMdl3LPV_IPCa');
OutTableON = array2table(OutDataLPV,'VariableNames',OutList);
[Won1,Kon1,refOn1] = getEvalCriteria(OutTableON,ref);

switchIPC = 2;
initDLRparameter(onIPC,switchIPC);
sim('test_SimulinkMdl3LPV_IPCa');
OutTableON = array2table(OutDataLPV,'VariableNames',OutList);
[Won2,Kon2,refOn2] = getEvalCriteria(OutTableON,ref);

switchIPC = 3;
initDLRparameter(onIPC,switchIPC);
sim('test_SimulinkMdl3LPV_IPCa');
OutTableON = array2table(OutDataLPV,'VariableNames',OutList);
[Won3,Kon3,refOn3] = getEvalCriteria(OutTableON,ref);

% IPC_K_1P  = 10^(-06);
% 
% IPC_KI_1p = 4.45e-05;
% IPC_KP_1p = -0.00010109;
% kp_Turm   = 1.9362e-05;
% IPC_KI_1p_tilt = IPC_KI_1p;
% IPC_KP_1p_tilt = IPC_KP_1p;


x0Ctrl = [IPC_K_1P, IPC_KI_1p,IPC_KP_1p,kp_Turm,IPC_KI_1p_tilt,IPC_KP_1p_tilt];
x0Ctrl  = zeros(1,6);

tic; [Won,Kon,refOn] = runSimIPC(x0Ctrl,ref,OutList); toc;

X = fminsearch(@(xCtrl) runSimIPC(xCtrl,ref,OutList),x0Ctrl);

save X X

options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
xSim = simulannealbnd(@(xCtrl) runSimIPC(xCtrl,ref,OutList),x0Ctrl,[],[],options);
