
close all; %clear; clc;
speedVec = [1,8,9,22]; %[1,2,8,9,21,22]; %[1,2,4,5,10,11,12,15,22]; % %18-3; %,9,22

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

cnt = 1;

aMatrix = nan(4*length(speedVec), 9);
for idxGapN = 1:length(fieldsGap)
    aFieldName = fieldsGap{idxGapN};
    aCellField = gapCell.(aFieldName);
    fprintf('\n %s\n---\n',aFieldName)
    for idx = speedVec
        str1 = sprintf('%2.2f & ',aCellField{idx});
        fprintf('%d m/s & %s \n',idx+3, str1);

        aMatrix(cnt,:) = aCellField{idx}(:);
        cnt = cnt +1;
    end
end


vec = 1:4;

figure;

for idx = 1 :9
    subplot(3,3,idx)
bar([aMatrix(vec,idx),aMatrix(vec+4,idx)])
end

figure;

for idx = 1 :9
    subplot(3,3,idx)
bar([aMatrix(vec+8,idx),aMatrix(vec+12,idx)])
end

%  nugap9DoF
% ---
% 4 m/s & 0.08 & 0.01 & 0.40 & 0.31 & 0.87 & 0.99 & 0.34 & 0.13 & 0.50 &  
% 11 m/s & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 &  
% 12 m/s & 0.03 & 0.03 & 0.36 & 0.99 & 1.00 & 1.00 & 0.35 & 0.23 & 0.83 &  
% 25 m/s & 0.01 & 0.08 & 0.23 & 0.96 & 0.85 & 1.00 & 0.06 & 0.23 & 0.77 &  
%  nugap9DoF
% ---
% 4 m/s & 0.08 & 0.01 & 0.40 & 0.31 & 0.87 & 0.99 & 0.34 & 0.13 & 0.50 &  
% 11 m/s & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 &  
% 12 m/s & 0.03 & 0.03 & 0.36 & 0.99 & 1.00 & 1.00 & 0.35 & 0.23 & 0.83 &  
% 25 m/s & 0.01 & 0.08 & 0.23 & 0.96 & 0.85 & 1.00 & 0.06 & 0.23 & 0.77 &  
% 
%  nugap5DoF
% ---
% 4 m/s & 0.10 & 0.01 & 0.41 & 0.31 & 0.97 & 0.99 & 0.34 & 0.43 & 0.56 &  
% 11 m/s & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 &  
% 12 m/s & 0.04 & 0.02 & 0.37 & 0.99 & 1.00 & 1.00 & 0.35 & 0.27 & 0.84 &  
% 25 m/s & 0.01 & 0.08 & 0.23 & 0.96 & 1.00 & 1.00 & 0.06 & 0.30 & 0.81 &  
% 
%  norm9DoF
% ---
% 4 m/s & 0.04 & 0.07 & 0.14 & 0.09 & 0.08 & 0.19 & 0.06 & 0.12 & 0.27 &  
% 11 m/s & 0.03 & 0.06 & 0.11 & 0.30 & 0.06 & 0.29 & 0.08 & 0.12 & 0.33 &  
% 12 m/s & 0.04 & 0.08 & 0.15 & 0.14 & 0.21 & 0.52 & 0.09 & 0.15 & 0.36 &  
% 25 m/s & 0.03 & 0.07 & 0.14 & 0.43 & 0.20 & 0.64 & 0.09 & 0.14 & 0.34 &  
%
%  norm9DoF
% ---
% 4 m/s & 0.04 & 0.07 & 0.14 & 0.09 & 0.08 & 0.19 & 0.06 & 0.12 & 0.27 &  
% 11 m/s & 0.03 & 0.06 & 0.11 & 0.30 & 0.06 & 0.29 & 0.08 & 0.12 & 0.33 &  
% 12 m/s & 0.04 & 0.08 & 0.15 & 0.14 & 0.21 & 0.52 & 0.09 & 0.15 & 0.36 &  
% 25 m/s & 0.03 & 0.07 & 0.14 & 0.43 & 0.20 & 0.64 & 0.09 & 0.14 & 0.34 &  
% % 
%  norm5DoF
% ---
% 4 m/s & 0.07 & 0.07 & 0.14 & 0.23 & 0.11 & 0.16 & 0.08 & 0.13 & 0.26 &  
% 11 m/s & 0.06 & 0.07 & 0.12 & 0.44 & 0.06 & 0.18 & 0.13 & 0.11 & 0.33 &  
% 12 m/s & 0.06 & 0.09 & 0.16 & 0.20 & 0.21 & 0.50 & 0.12 & 0.15 & 0.36 &  
% 25 m/s & 0.04 & 0.06 & 0.18 & 0.50 & 0.23 & 0.60 & 0.13 & 0.17 & 0.32 &  
