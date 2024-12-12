
close all; %clear; clc;
speedVec = [1,8,9,22]; %[1,2,8,9,21,22]; %[1,2,4,5,10,11,12,15,22]; % %18-3; %,9,22

useActuatorStates = 0;
figNoAdd = 0;
createBodePlots = 1;
plotVisible = 1;
noOut = 3;

figFolder = 'figLinDiss02';

%compareLinearModels(speedVec,figFolder, useActuatorStates,figNoAdd,createBodePlots,plotVisible)
[sysOut,gapCell] = compareLinearModels2(speedVec,figFolder,useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut);

%% Read out the quantative information

fieldsGapAll = fieldnames(gapCell); % get the field names
fieldsGap = fieldsGapAll(contains(fieldsGapAll,'n')); % only use Vinnicombe/nu-gap and norm for now
sortFieldsGap = sort(fieldsGap); % Mdl1 (5 states) should be listed before Mdl2 (9 states)

cnt = 1; % this is the counter for constructing
aMatrix = nan(4*length(speedVec), 9); %the matrix with the values

for idxGapN = 1:length(fieldsGap)
    aFieldName = sortFieldsGap{idxGapN};
    aCellField = gapCell.(aFieldName);
    fprintf('\n %s\n---\n',aFieldName)

    for idx = speedVec
        str1 = sprintf('%2.2f & ',aCellField{idx});
        fprintf('%d m/s & %s \n',idx+3, str1); % simply write the tables out
        tmp = aCellField{idx}'; %Transpose: Matlab is column first
        aMatrix(cnt,:) = tmp(:); % but the subplots are row first
        cnt = cnt + 1; % next row of the matrix
    end
end


%% Set informations for both plots

lc = lines; %line colors
lc1 = [0*[0.5, 0.5,0.5];lc];
len = length(speedVec);
vec = 1:len;
idxVec = 1:9; % this is the number of input/output combination

tmpStr = sprintf('%d m/s,', speedVec + 3); %wind speed infor
windTickLabel = regexp(tmpStr(1:end-1), ',','split');

% sysOutputname =[{'RotSpd \omega_r (rad/s)'};{'TwrAcc_{fa} (m/s^2)'}; {'TwrAcc_{sw} (m/s^2)'}];
% sysInputname =[{'GenTq T_g (kNm)'}; {'BlPitch \beta_0 (rad)'};{'Wind V_{\infty} (m/s)'}];

sysOutputname =[{'RotSpd \omega_r (-)'};{'TwrAcc_{fa} (-)'}; {'TwrAcc_{sw} (-)'}];
sysInputname =[{'GenTq T_g (-)'}; {'BlPitch \beta_0 (-)'};{'Wind V_{\infty} (-)'}];

Mdl2Str = ' Rot,Twr,Gen+Bld';


%% Set values for first, norm plot
aMat1 = [aMatrix(vec,idxVec),aMatrix(vec+len,idxVec)];
aMax = max(aMat1(:))* 1.01;
legCell = ['                        Norm: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

figure(1);
for idx = 1:9
    subplot(3,3,idx)
    b = bar([aMatrix(vec,idx),aMatrix(vec+len,idx)]);

    for k = 1:length(b), b(k).FaceColor = lc1(k,:); end
    set(gca,'XTickLabel',windTickLabel);
    set(gca,'YLim',[0,aMax]);

    if idx == 1
        title( legCell)
    end

    if mod(idx-1,3)== 0
        anOutputStr = sysOutputname{(idx-1)/3 + 1};
        ylabel(anOutputStr);
    end

    if idx >= 7
        anInputStr = sysInputname{idx-6};
        xlabel(anInputStr);
    end
end

figStr = 'Norm';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');


%% Set values for second, nugap plot
aMat1 = [aMatrix(vec+2*len,idxVec),aMatrix(vec+3*len,idxVec)];
aMax = max(aMat1(:))* 1.01;
legCell = ['                          Nugap: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

figure(2);
for idx = 1:9
    subplot(3,3,idx)
    b = bar([aMatrix(vec+ 2*len,idx),aMatrix(vec+3*len,idx)]);

    for k = 1:length(b), b(k).FaceColor = lc1(k,:); end
    set(gca,'XTickLabel',windTickLabel);
    set(gca,'YLim',[0,aMax]);

    if idx == 1
        title( legCell)
    end

    if mod(idx-1,3)== 0
        anOutputStr = sysOutputname{(idx-1)/3 + 1};
        ylabel(anOutputStr);
    end

    if idx >= 7
        anInputStr = sysInputname{idx-6};
        xlabel(anInputStr);
    end
end

figStr = 'Nugap';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');


%%
aTest = reshape(1:9,3,3);
aTest = aTest';
aTestVec = aTest(:);

aPerm5DOF = permute(gapCell.mag5DoF{1},[2 1 3]); % abs(permute(gapCell.mag5DoF{1},[2 1 3]) - permute(gapCell.magFAST{1},[2 1 3]));
aPerm9DOF = permute(gapCell.mag9DoF{1},[2 1 3]); %abs(permute(gapCell.mag9DoF{1},[2 1 3]) - permute(gapCell.magFAST{1},[2 1 3]));

figure;

for idx = 1:9

    idxOut = mod((idx-1),3) + 1;
    idxIn = floor((idx-1)/3) + 1;

    mag_9DOF = 20*log10(squeeze(aPerm9DOF(idxOut,idxIn,:)));    % Squeeze to get a 1D array
    mag_5DOF = 20*log10(squeeze(aPerm5DOF(idxOut,idxIn,:)));
    % mag_FAST = 20*log10(squeeze(magFAST(idxOut,idxIn,:)));

    subplot(3,3,idx)
    semilogx(gapCell.freq9DoF{1},mag_9DOF,'-',...
        gapCell.freq9DoF{1},mag_5DOF,'k--');
    axis tight, grid on;

    if mod(idx-1,3)== 0
        anOutputStr = sysOutputname{(idx-1)/3 + 1};
        ylabel(anOutputStr);
    end

    if idx >= 7
        anInputStr = sysInputname{idx-6};
        xlabel(anInputStr);
    end

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
%  norm9DoF
% ---
% 4 m/s & 0.05 & 0.07 & 0.14 & 0.09 & 0.08 & 0.19 & 0.06 & 0.12 & 0.27 &
% 11 m/s & 0.04 & 0.06 & 0.11 & 0.30 & 0.06 & 0.29 & 0.08 & 0.12 & 0.33 &
% 12 m/s & 0.04 & 0.08 & 0.15 & 0.14 & 0.21 & 0.52 & 0.09 & 0.15 & 0.36 &
% 25 m/s & 0.04 & 0.06 & 0.17 & 0.43 & 0.20 & 0.64 & 0.09 & 0.14 & 0.34 &

%  norm5DoF
% ---
% 4 m/s & 0.07 & 0.07 & 0.14 & 0.23 & 0.11 & 0.16 & 0.08 & 0.13 & 0.26 &
% 11 m/s & 0.06 & 0.07 & 0.12 & 0.44 & 0.06 & 0.18 & 0.13 & 0.11 & 0.33 &
% 12 m/s & 0.06 & 0.09 & 0.16 & 0.20 & 0.21 & 0.50 & 0.12 & 0.15 & 0.36 &
% 25 m/s & 0.04 & 0.06 & 0.18 & 0.50 & 0.23 & 0.60 & 0.13 & 0.17 & 0.32 &
