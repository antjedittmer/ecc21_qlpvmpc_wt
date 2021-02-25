clc; clear; bdclose;

% Add Paths for simulation
addpath('MdlWECS');
addpath('subfunctions');
addpath('subfunctionsMPC');

% initialize constants to base workspace for CPC simulation
testCtrl;
load(fullfile(dataInPath,'OutDataWind18NTM.mat'),'OutData','OutTable');

onIPC = 0;
switchIPC = 1;
initDLRparameter(onIPC,switchIPC);
[sysd,x0,uoff,xoff,yoff] = initLPVParam;

%% Get the baseline values with PI CPC
if exist('OutTableWeightsOrig.mat','file') == 2
    load('OutTableWeightsOrig.mat','W','K','ref','weightVec')
else
    sim('test_SimulinkMdl3LPV_IPCa');
    OutTableCPC = array2table(OutDataLPV,'VariableNames',OutList);
    
    weightVec = [1.25,1,0.5,1,1.25];
    [W,K,ref,weightVec] = getEvalCriteria(OutTableCPC,[],weightVec);
    save('OutTableWeightsOrig.mat','W','K','ref','weightVec')
end

% plotTimeseriesCtrlComparison(idxTime, outTableNRELorig, outTableNRELnew,outTableDLR, strNewCtrl,figNo1,figNo2,testCaseStr)
% plotTimeseriesCtrlComparison([], OutTableCPC, OutTableCPC,OutTableCPC)

onIPC = 1;
switchIPC = 1;
initDLRparameter(onIPC,switchIPC);

% IPC optimization with X
if exist('X.mat','file') == 2
    load('X','X')
else
    x0Ctrl  = zeros(1,3);
    X = fminsearch(@(xCtrl) runSimIPC(x0Ctrl,ref,OutList),x0Ctrl);
    save('X','X');
end
[WIPC,KIPC,refOPC,OutTableIPC] = runSimIPC(X,ref,OutList,weightVec);

if exist('xSim.mat','file') == 2
    load('xSim','xSim');
else
    options = optimoptions('simulannealbnd','PlotFcns',...
        {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
    xSim = simulannealbnd(@(xCtrl) runSimIPC(xCtrl,ref,OutList),x0Ctrl,[],[],options);
    save('xSim','xSim');
end

[WIPCSimAn,KIPCSimAn,refOPCSimAn,OutTableIPCSimAn] = runSimIPC(xSim,ref,OutList,weightVec);

weightVec = [1,1,1,1,100];
[W,K,ref100Tower,weightVec100Tower] = getEvalCriteria(OutTableCPC,[],weightVec);
 X100 = fminsearch(@(xCtrl) runSimIPC(x0Ctrl,ref100Tower,OutList,weightVec100Tower),x0Ctrl);
 
strNewCtrl = 'IPC_{SimAnn}';
plotTimeseriesCtrlComparison([], OutTableCPC,OutTableIPC,OutTableIPCSimAn,strNewCtrl)


cellOfTable{1} = OutTableCPC;
cellOfTable{2} = OutTableIPC;
cellOfTable{3} = OutTableIPCSimAn;
calculateEvalCriteria(cellOfTable)


varnames = {'Wind', 'RotSpeed', 'GenPwr', 'BlPitch1', 'GenTq', ...
    'NcIMUTAxs', 'NcIMUTAys'};

OutTableMPC = array2table(OutDataTest_4to25,'VariableNames',varnames);
strNewCtrl = 'qLPV MPC';
plotTimeseriesCtrlComparison([], OutTableCPC,OutTableIPC,OutTableMPC,strNewCtrl)


%% Unused code; Check results from thesis
% onIPC = 1;
% switchIPC = 1;
% initDLRparameter(onIPC,switchIPC);
% sim('test_SimulinkMdl3LPV_IPCa');
% OutTableON = array2table(OutDataLPV,'VariableNames',OutList);
% [Won1,Kon1,refOn1] = getEvalCriteria(OutTableON,ref);
% 
% switchIPC = 2;
% initDLRparameter(onIPC,switchIPC);
% sim('test_SimulinkMdl3LPV_IPCa');
% OutTableON = array2table(OutDataLPV,'VariableNames',OutList);
% [Won2,Kon2,refOn2] = getEvalCriteria(OutTableON,ref);
% 
% switchIPC = 3;
% initDLRparameter(onIPC,switchIPC);
% sim('test_SimulinkMdl3LPV_IPCa');
% OutTableON = array2table(OutDataLPV,'VariableNames',OutList);
% [Won3,Kon3,refOn3] = getEvalCriteria(OutTableON,ref);
% 
