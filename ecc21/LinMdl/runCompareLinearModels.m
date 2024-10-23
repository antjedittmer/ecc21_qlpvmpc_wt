
close all; clear; clc;
speedVec = [1,8,9,22]; %[1,2,8,9,21,22]; %[1,2,4,5,10,11,12,15,22]; % %18-3;

useActuatorStates = 0;
figNoAdd = 0;
createBodePlots = 1;
plotVisible = 1;
noOut = 3;

figFolder = 'figLinUIO1';

%compareLinearModels(speedVec,figFolder, useActuatorStates,figNoAdd,createBodePlots,plotVisible)
[sysOut,gapCell] = compareLinearModels(speedVec,figFolder,useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut);

fieldsGapAll = fieldnames(gapCell);
fieldsGap = fieldsGapAll(contains(fieldsGapAll,'n')); % only use Vinnicombe gap and norm for now

for idxGapN = 1:length(fieldsGap)
    aFieldName = fieldsGap{idxGapN};
    aCellField = gapCell.(aFieldName);
    fprintf('\n %s\n---\n',aFieldName)
    for idx = speedVec
        str1 = sprintf('%2.2f & ',aCellField{idx});
        fprintf('%d m/s & %s \n',idx, str1);
    end
end