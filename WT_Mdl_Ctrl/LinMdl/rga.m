clear; clc; close all;

%% Initialize name of models and load into workspace
addpath(fullfile(fileparts(pwd),'dataIn'));
idxModelNames = 1; %1:5
modelNames =[{'NREL5MW_linearised_Gen_Fa_SS_Flap1'};
    {'NREL5MW_linearised_Gen_Fa_Flap1'};
    {'NREL5MW_linearised'};
    {'NREL5MW_linearised_Gen_DriveTrain_Fa_SS_Flap1'}; {'NREL5MW_linearised'} ];
modelNames = modelNames(idxModelNames);

modelsL = [{'Fast Rot+FA+SW+1^{st} flap mode (12 states)'};
    {'FAST 10 states'};
    {'FAST 30 states'};
    {'Fast Rot+drive train+FA+SW+1^{st} flap mode (14 states)'};
    {'Fast Rot+FA (4 states)'}];

modelsL = modelsL(idxModelNames);

lenModelNames = length(idxModelNames);
sysmCell = cell(lenModelNames,1);
sysmCell_MBC = cell(lenModelNames,1);
for idxModels = 1: lenModelNames
    load(modelNames{idxModels},'sysm','sysm_MBC');
    sysmCell{idxModels} = sysm;
    sysmCell_MBC{idxModels} = sysm_MBC;
end

lowerb = -2;
upperb = 1.5;
V = 4:25;
for index = 1:21 %[1:4,9,10,21]
    for idxModelNames=1:lenModelNames
        sysm = sysmCell{idxModelNames};
        % model = [1 0 0; 0 1 0; 0 0 2*pi/60]*sysm{index}([12 6 38],[8 9 1])*[1000 0 0; 0 1 0; 0 0 1]; %pi/180
        model = [2*pi/60 0; 0 1000/97] * sysm{index}([33 39],[9 8])*[1 0; 0 1]; %pi/180
        model.OutputName = sysm{index}.OutputName([33 39]);
        model.InputName = sysm{index}.InputName([9 8]);
        %bodemag(model,opts); hold on
        G = tf(model);
        titlestr = sprintf(' Wind speed %d m/s', V(index));
        lambda = plot_rga(G,lowerb,upperb,titlestr);
        grid on
    end
end
