clc; clear; bdclose;

% Add Paths for simulation
addpath('WECSMdl');
addpath('subfunctions');
addpath('subfunctionsMPC');
addpath('subfunctionsUtil');
% initialize constants to base workspace for CPC simulation
testCtrl;
load(fullfile(dataInPath,'OutDataWind18NTM.mat'),'OutData','OutTable');

outDataDir = 'dataOut3';

if ~isempty(outDataDir) && ~isdir(outDataDir) %#ok<ISDIR>
    mkdir(outDataDir);
end

loadMat = 1;
onIPC = 0;
switchIPC = 1;
initDLRparameter(onIPC,switchIPC);
[sysd,x0,uoff,xoff,yoff] = initLPVParam;

%% Get the baseline values with PI CPC
matWeightsOrig = fullfile(outDataDir,'OutTableWeightsOrig.mat');
if loadMat && exist(matWeightsOrig,'file') == 2
    load(matWeightsOrig,'W','K','ref','weightVec','OutTableCPC')
else
    sim('test_SimulinkMdl3LPV_IPCa');
    OutTableCPC = array2table(OutDataLPV,'VariableNames',OutList);
    
    weightVec = [50,1,0.5,5,5];
    [W,K,ref,weightVec] = getEvalCriteria(OutTableCPC,[],weightVec);
    save(matWeightsOrig,'W','K','ref','weightVec','OutTableCPC')
end

% plotTimeseriesCtrlComparison(idxTime, outTableNRELorig, outTableNRELnew,outTableDLR, strNewCtrl,figNo1,figNo2,testCaseStr)
% plotTimeseriesCtrlComparison([], OutTableCPC, OutTableCPC,OutTableCPC)

onIPC = 1;
switchIPC = 1;
initDLRparameter(onIPC,switchIPC);

% IPC optimization with X
matWeightsIPCfminsearch = fullfile(outDataDir,'X.mat');
if loadMat && exist(matWeightsIPCfminsearch,'file') == 2
    load(matWeightsIPCfminsearch,'X')
else
    x0Ctrl  = zeros(1,3);
    X = fminsearch(@(xCtrl) runSimIPC(xCtrl,ref,OutList,weightVec),x0Ctrl);
    save(matWeightsIPCfminsearch,'X');
end
[WIPC,KIPC,refOPC,OutTableIPC] = runSimIPC(X,ref,OutList,weightVec);


x = ga( runSimIPC,nvars,A,b,[],[],lb,ub,nonlcon,IntCon,options)

%% Set vector which pitch rate constraints to be test
maxRateVector = [13,10,8];

%% Loop over pitch rate constraint vector:
nRate = length(maxRateVector);
currentOutTableTestCell = cell(nRate,1); % cell with

%dataDirOut = 'dataOut';

DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);

currentOutTableIPCCell = cell(nRate,1); % cell with

simMdlname = 'test_SimulinkMdl2_Baseline';
x0Ctrl = X;
for idx = 1: nRate
    maxRate = maxRateVector(idx);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlCtrlObj,controlValue);
    saveChanges(myDictionaryCtrlObj)
    
    matFileOutTableTest1 = fullfile(outDataDir,...
        sprintf('OutTableTest_rate%02d_IPC.mat',maxRate));
    
    
    if exist(matFileOutTableTest1,'file') == 2
        load(matFileOutTableTest1,'xSim');
    else
        xSim = simulannealbnd(@(xCtrl) runSimIPC(xCtrl,ref,OutList,weightVec),x0Ctrl);
        save(matFileOutTableTest1,'xSim');
    end
    currentOutTableIPCCell{idx} = xSim;
end
load('xSim','xSim');

[WIPCSimAn,KIPCSimAn,refOPCSimAn,OutTableIPCSimAn] = runSimIPC(xSim,ref,OutList,weightVec);

%weightVec = [1,1,1,1,100];
[W,K,ref100Tower,weightVec100Tower] = getEvalCriteria(OutTableCPC,[],weightVec);
[WIPC,KIPC,ref100TowerIPC,weightVec100TowerIPC] = getEvalCriteria(OutTableIPC,[],weightVec);
x0Ctrl  = zeros(1,3);
%X100 = fminsearch(@(xCtrl) runSimIPC(x0Ctrl,ref100Tower,OutList,weightVec100Tower),x0Ctrl);

strNewCtrl = 'IPC_{SimAnn}';
plotTimeseriesCtrlComparison([],OutTableCPC,OutTableIPC,OutTableIPCSimAn,strNewCtrl)

%[K_matrix, strK] = calculateEvalCriteria(aTestCell3,summaryTickLabel);


cellOfTable{1} = OutTableCPC;
cellOfTable{2} = OutTableIPC;
cellOfTable{3} = OutTableIPCSimAn;
[K_matrix, strK] = calculateEvalCriteria(cellOfTable);

varnames = {'Wind', 'RotSpeed', 'GenPwr', 'BlPitch1', 'GenTq', ...
    'NcIMUTAxs', 'NcIMUTAys'};


load('OutTableTest_rate13_MPC_LPV.mat','OutTableTest_rate13_MPC');
load('OutTableTest_rate01_MPC_LPV.mat','OutTableTest_rate01_MPC');
load('OutTableTest_rate04_MPC_LPV.mat','OutTableTest_rate04_MPC');
load('OutTableTest_rate02_MPC_LPV.mat','OutTableTest_rate02_MPC');
load('OutTableTest_rate08_MPC_LPV.mat','OutTableTest_rate08_MPC');
% OutTableMPC = OutTableTest_rate13_MPC; %array2table(OutDataTest_4to25,'VariableNames',varnames);
% OutTableMPC.zeta_dotdot = zeta_dotdot;
% strNewCtrl = 'qLPV MPC';
% 
% summaryTickLabel = {'CPC PI', 'MPC'};
% 
% OutTablePICPC = array2table(OutDataTest,'VariableNames',[varnames,'zeta_dotdot']);
% OutTablePICPC.Time = OutTableCPC.Time;
cellOfTableNew{1} = cellOfTable{1};
cellOfTableNew{2} = OutTableIPC;
cellOfTableNew{3} = OutTableTest_rate13_MPC;
% cellOfTableNew{4} = OutTableTest_rate08_MPC;
% cellOfTableNew{5} = OutTableTest_rate04_MPC;
cellOfTableNew{4} = OutTableTest_rate01_MPC;

summaryTickLabel = {'PI, limit 8 deg/s', 'IPC, limit 8 deg/s','MPC, limit 13 deg/s', 'MPC, limit 1 deg/s'};
calculateEvalCriteria(cellOfTableNew,summaryTickLabel)

plotTimeseriesCtrlComparison([], OutTableCPC,OutTableIPC,OutTableIPC,strNewCtrl)


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
