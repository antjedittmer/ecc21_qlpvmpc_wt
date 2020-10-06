%% Arm-Driven Inverted Pendulum simulation model simADIP
% Display simulation model and results used as basis for the paper 
% Pablo S.G. Cisneros, Herbert Werner, 'Wide Range Stabilization of a 
% Pendubot using quasi-LPV Predictive Control', 3rd IFAC Workshop on Linear 
% Parameter Varying Systems LPVS 2019: Eindhoven, Netherlands,4–6 Nov. 2019

clc; clear; close all; bdclose all;
loadSimResults = 0;

%% Get Simulation model name
% Get the name of the model. Model is saved with a different name in R2018b.
simFile = '.slx';
simulinkModel = ['simADIP',simFile]; 
versionRelease = version('-release');

is2018b = strcmp('2018b',versionRelease);
if is2018b 
    simulinkModel = strrep(simulinkModel, simFile,['_',versionRelease,simFile]);
end

% Set outputfolder names
outDirFig = ['outFig',versionRelease];
if ~isdir(outDirFig) %#ok<ISDIR> isdir used for R2007b 
    mkdir(outDirFig)
end

outDirMat = ['outMat',versionRelease];
if ~isdir(outDirMat) %#ok<ISDIR> isdir used for R2007b 
    mkdir(outDirMat)
end
simoutFile = fullfile(outDirMat,'simout.mat');

%% Get scope names
% Get handle on scopes and set their data logging on

% Load Simulink model
simulinkModelPrint = simulinkModel(1:end-4);
open_system(simulinkModelPrint);

% Get handle on all scopes and print out
aBlockTypeScope = 'Scope';
listScopes = find_system(simulinkModelPrint,...
    'LookUnderMasks','all','BlockType', aBlockTypeScope);
fprintf('%s\n',listScopes{:});


%% Run Simulation model or load simulation data
% Run Simulation model or load simulation data

% Run Simulation model if data is not available
if ~(exist(simoutFile,'file')== 2) || ~loadSimResults   
    % Set logging on for all scopes
    for idxScopes = 1: length(listScopes)
        aScope = listScopes{idxScopes};
        set_param(aScope,'SaveToWorkspace','on', 'DataFormat','StructureWithTime',...
        'DataLoggingVariableName',['ScopeData',num2str(idxScopes)])
    end
    sim(simulinkModel);
    save(simoutFile,'ScopeData*')
    
    % Clean up
    clear ScopeData* 
    for idxScopes = 1: length(listScopes)
        aScope = listScopes{idxScopes};
        set_param(aScope,'SaveToWorkspace','off')
    end
end
simData = load(simoutFile);


%% Display the Simulink model as image
% Print Simulink model to png file. Color inports and outports to make them
% easily detectable.

% Get name of model, color in/ouports, create outDirfolder
colorInOutports(simulinkModelPrint ,'cyan','yellow');

% create full path to png file and display in MATLAB
simulinModelFile = fullfile(outDirFig,[simulinkModelPrint,'.png']);
print(['-s',simulinkModelPrint],'-dpng',simulinModelFile);
img = imread(simulinModelFile);

if is2018b
    figure(1);
    imshow(img);
end

%% Display figures main level
% Display scopes of main level in one plot.

% Gets the scopes on base layer 
simDataCell = struct2cell(simData);
simDataName = cellfun(@(x) x.blockName, simDataCell, 'UniformOutput',false);
cellScopesSim = regexp(simDataName, '/', 'split');
isMainScope = cellfun(@length,cellScopesSim) == 2; 

simDataCellMain = simDataCell(isMainScope);
simDataStructSub = simDataCell(~isMainScope);
simDataCellSubName = cellScopesSim(~isMainScope);

lenBlockName = cell2mat(cellfun(@(x) max(strfind(x,'/')),...
    simDataName(~isMainScope), 'UniformOutput',false));

% create the plots
% hf = figure(is2018b + 1); % there is a figure in 2018b, no figure in 2008b
for idxMain = 1: length(simDataCellMain)
    aScopeStruct = simDataCellMain{idxMain}; %  list is ordered alphatically
    %subplot(length(simDataCellMain),1,idxMain)
    hf = figure(is2018b + idxMain);
    nrSubplots = length(aScopeStruct.signals);
    
    for idxSub = 1: nrSubplots
        subplot( nrSubplots,1,idxSub)
        plot(aScopeStruct.time,aScopeStruct.signals(idxSub).values);
        axis tight; grid on;
        
        if idxSub == 1
            xlabel('Time (s)');
            ht = title(aScopeStruct.blockName); %ht.Interpreter = 'none';
            set(ht,'Interpreter','none');
        end
    end
    
    
    if isnumeric(hf)
        saveas(hf,fullfile(outDirFig,sprintf('Scope%02d',hf)),'png');
    end
end

%% Get subsystem names
% Get subsystem names with find_system. Print out with fprintf.
% Remove subsystems without blocks from list.
%
% * Simulation/MATLAB Function'
% * Simulation/MATLAB Function2'
% * Simulation/Signal Builder'

% Find subsystems with find_system
aBlockType = 'SubSystem';
listSubsysAll = find_system(simulinkModelPrint,...
    'LookUnderMasks','all','BlockType', aBlockType);

% Get listSubsys for user-build subsystems
listSubsys = setdiff(listSubsysAll,...
    {[simulinkModelPrint, '/MATLAB Function'],...
    [simulinkModelPrint, '/MATLAB Function2'],...
    [simulinkModelPrint, '/Signal Builder']});

fprintf('%s\n',listSubsys{:});


%% Display subsystems as images
% Displaye subsystems as images with imshow (not available in R2007b)
% If there are only a  be few subsystems displaying them by looping
% over the list listSubsys is a valid, easily readable solution. 

for idxSubsys = 1: length(listSubsys)
    aSimulinkModelPrint = listSubsys{idxSubsys};
    fprintf('%s\n',aSimulinkModelPrint);
    if ~is2018b
        open_system(aSimulinkModelPrint); 
    end
    aPic = fullfile(outDirFig,sprintf('%s_Subsys%02d.png',simulinkModelPrint,idxSubsys));
    print(['-s', aSimulinkModelPrint],'-dpng',aPic)
    if ~is2018b
        close_system(aSimulinkModelPrint);
    end
    
    if is2018b
        img = imread(aPic);
        stringTitle = strrep(aSimulinkModelPrint,[simulinkModelPrint,'/'],'');
        hf = figure(hf.Number+1);
        imshow(img);
        ht = title(stringTitle); %ht.Interpreter ='none';
        set(ht,'Interpreter','none');
    end
     
    drawnow;
end

%% Display content of MATLAB functions
% Display content of MATLAB functions: Get sfroot.

% Get handle to MATLAB functions of simulation model
S = sfroot;
stateflowEMCharts = S.find('-isa','Stateflow.EMChart');

% Display block name and content
for idxSF = 1: length(stateflowEMCharts)
    aMatlabFcn = stateflowEMCharts(idxSF);
    disp(aMatlabFcn.Path)
    fprintf('---\n')
    disp(aMatlabFcn.Script)
    fprintf('\n\n')
end


