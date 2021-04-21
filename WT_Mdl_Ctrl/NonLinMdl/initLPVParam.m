function [sysd,x0,uoff,xoff,yoff] = initLPVParam(DT,tmpLinStr,OutList,OutTable)
%addpath(fullfile(fileparts(pwd),'dataIn'));
% initModel5MWNREL;
% load('workSpace.mat')
% load('OutDataSweep.mat');
% OutTable = array2table(OutData,'VariableNames', OutList);
% idxWind = contains(OutList, 'Wind');
% vectWind = OutTable{:,idxWind};
% vectAmpWind = sqrt(sum((vectWind.^2),2));
if ~nargin
    DT = 0.0080;
end

if nargin < 2 || isempty(tmpLinStr)
    tmpLinStr = fullfile(fileparts(pwd),'dataIn','NREL5MW_linearised_4to25.mat');
end

if nargin < 3 || isempty(OutList)
    load('workSpace.mat','OutList');
end

if nargin < 4 || isempty(OutList)
    load('workSpace.mat','OutList');
end

if nargin < 5 || isempty(OutTable)
    % Orig code: load('OutDataSweep.mat','OutData'); OutTable = array2table(OutData,'VariableNames', OutList);
    load('OutDataSweep.mat','OutTable'); %OutTable = array2table(OutData,'VariableNames', OutList);
end
v0 = OutTable.Wind1VelX(1);

% Load linearized file
tmpLin = load(tmpLinStr);

% Order output, select input
ssTestOutputName = tmpLin.sysm{1}.OutputName;
filterInput = [1,0,0,1,1,1,0,1,0] == 1;
filterOutput = NaN(length(ssTestOutputName),1);
for idxOut = 2 : length(OutList)
    noOut = find(contains(tmpLin.sysm{1}.OutputName,OutList(idxOut)));
    filterOutput(idxOut-1) = noOut;
end

% Create state space array
sizeSS = size(tmpLin.sysm{1});
sizeSS1= length(tmpLin.sysm{1}.A);

sysd = drss(sizeSS1,sizeSS(1),sum(filterInput),length(tmpLin.sysm),1);
for idx = 1: length(tmpLin.sysm)
    ssTest = tmpLin.sysm{idx};

    ssTestIn = c2d(ss(ssTest.A,ssTest.B(:,filterInput),ssTest.C(filterOutput,:),...
        ssTest.D(filterOutput,filterInput)),...
        DT,'zoh');
    ssTestIn.InputName = ssTest.InputName(filterInput);
    ssTestIn.OutputName = ssTest.OutputName;
    ssTestIn.StateName = ssTest.StateName;
    sysd(:,:,idx,1) = ssTestIn;
end
sysd.SamplingGrid =  struct('V',tmpLin.Lin.V);

% Set equilibrium points
uoff1 = tmpLin.Lin.u_op; %cellfun(@(x) x(x>0), C, 'UniformOutput',false);
uoff2 = cell2mat(cellfun(@(x) x(filterInput), uoff1,'UniformOutput', false));
uoff = reshape(uoff2,size(uoff2,1),1, size(uoff2,2)); 

xoff2 = cell2mat(tmpLin.Lin.x_op);
xoff = reshape(xoff2,size(xoff2,1),1, size(xoff2,2));

yoff1 = tmpLin.Lin.y_op;
yoff2 = cell2mat(cellfun(@(x) x(filterOutput), yoff1,'UniformOutput', false));
yoff = reshape(yoff2,size(yoff2,1),1, size(yoff2,2));

x0 = xoff(:,min(abs(tmpLin.Lin.V -  v0)) == abs(tmpLin.Lin.V -  v0));
