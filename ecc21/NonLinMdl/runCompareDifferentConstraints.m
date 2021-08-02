clc;
initWorkspace;


%% Set path and parameter names
simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta'; 
outDataSimulationMat = 'OutDataWind18NTW.mat';
loadData = 1;

% Set path to data output directory
workDir = fileparts(mfilename('fullpath'));
dataDirOut = fullfile(workDir, 'dataOut');
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

%% Get Data Directory value of pitch actuator rate
DDName = 'DD_MdlCtrl_qLPVMPC.sldd';
myDictionaryObj = Simulink.data.dictionary.open(DDName);
dDataSectObj = getSection(myDictionaryObj,'Design Data');
controlObj = getEntry(dDataSectObj,'Control');
controlValue = getValue(controlObj);

%% Set vector which pitch rate constraints to be test
maxRateVector = 1:13;

%% Loop over pitch rate constraint vector: 
nRate = length(maxRateVector);
currentOutTableTestCell = cell(nRate,1); % cell with
tictoc_LPVMPCcell = cell(nRate,1);
GenPwrRefCell= cell(nRate,1);

simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta'; 


for idx = 1: nRate
    maxRate = maxRateVector(idx);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlObj,controlValue);
    saveChanges(myDictionaryObj)
    
    matFileOutTableTest1 = fullfile(dataDirOut,...
        sprintf('OutTableTest_rate%02d_MPC.mat',maxRate));
    
    [currentOutTableTestCell{idx}, tictoc_LPVMPCcell{idx},GenPwrRefCell{idx}] = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname);
end

% Plot qLmpc result
OutTable4 = currentOutTableTestCell{1};
OutTable4.Time = OutTable.Time(1:height(OutTable4));
testCaseStr =  'qLMPC CPC Rate Constraints [deg/s]';
testCaseCell = {'standard (8) ','high (13) ', 'low (4) ',};
plotOutTable(OutTable4,currentOutTableTestCell{5},currentOutTableTestCell{end},...
   testCaseStr,testCaseCell);


%% Get Data Directory value of pitch actuator rate
DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);

currentOutTablePICell = cell(nRate,1); % cell with

simMdlname = 'test_SimulinkMdl2_Baseline';

for idx = 1: nRate
    maxRate = maxRateVector(idx);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlCtrlObj,controlValue);
    saveChanges(myDictionaryCtrlObj)
        
    matFileOutTableTest1 = fullfile(dataDirOut,...
        sprintf('OutTableTest_rate%02d_PI.mat',maxRate));
    
    currentOutTablePICell{idx} = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname);
end

% Plot qLmpc result
OutTable4 = currentOutTablePICell{1};
OutTable4.Time = OutTable.Time(1:height(OutTable4));
testCaseStr =  'PI CPC Rate Constraints [deg/s]';
figNo2 = 2;
figDir = '';
strFig = 'testConstrPI';
plotOutTable(OutTable4,currentOutTablePICell{5},currentOutTablePICell{end},...
   testCaseStr,testCaseCell, figNo2,figDir,strFig)

%% Extract signals from table
idxT = 60/DT:height(currentOutTablePICell{1});
matrixPitchPI = cell2mat(cellfun(@(x) x.BlPitch1(idxT)', currentOutTablePICell, 'UniformOutput',false));
matrixPitchQLMPC = cell2mat(cellfun(@(x) x.BlPitch1(idxT)', currentOutTableTestCell, 'UniformOutput',false));
matrixTwrPI = cell2mat(cellfun(@(x) x.NcIMUTAxs(idxT)', currentOutTablePICell, 'UniformOutput',false));
matrixTwrQLMPC = cell2mat(cellfun(@(x) x.NcIMUTAxs(idxT)', currentOutTableTestCell, 'UniformOutput',false));

diffPitchPI = diff(matrixPitchPI');
diffPitchQLMPC = diff(matrixPitchQLMPC');

stdDiffPitchPI = std(abs(diff(matrixPitchPI')))*180/pi;
stdDiffPitchQLMPC = std(abs(diff(matrixPitchqLMPC')))*180/pi;
stdTowerPI = std(abs(matrixTwrPI'));
stdTowerQLMPC = std(abs(matrixTwrQLMPC'));

figure; 
subplot(2,1,1);
plot(maxRateVector,stdDiffPitchPI,'*-',...
    maxRateVector,stdDiffPitchQLMPC,'o-');
legend('PI','qLMPC');
axis tight; grid on;
ylabel('std_{\Delta Pitch} [deg/s]');
subplot(2,1,2);
plot(maxRateVector,stdTowerPI,'*-',...
    maxRateVector,stdTowerQLMPC,'o-');
axis tight; grid on;
ylabel('std_{\Delta Pitch} [deg/s]');
xlabel('std_{TwrFA} [m/s^2]');
figure;
plot(stdDiffPitchPI,stdTowerPI,'*-',...
stdDiffPitchQLMPC,stdTowerQLMPC,'o-');



