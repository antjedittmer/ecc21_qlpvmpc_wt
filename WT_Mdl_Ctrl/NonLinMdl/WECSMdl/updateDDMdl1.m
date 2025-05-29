% script to update data dictionary
% clear;
function Lin_points = updateDDMdl1(multRb,useTurn, x)

if ~nargin
    multRb = 0.75;
end

if  nargin < 2
    useTurn = 0;
end

if nargin <3
    x = [1,1];
end


% Update DD
nameCell = {'beta';'Ce';'Cq';'Ct';'K5'; 'lambda'; 'M';'Q'; 'rho';'wecs';'Lin_points'};

loadEqPoints = 1; % load equilibrium/reference points for wind speeds from file
% if  loadEqPoints == 1
dataInLin = fullfile(pwd,'dataIn','Lin_points1.mat');
load(dataInLin,'Lin'); %load('NREL5MW_linearised_4to25.mat','Lin');

%else %get them from steady state of long simulation run
load('OutDataStep.mat', 'OutTable');
idxT = OutTable.Time >= 100;
OutTableT = OutTable(idxT,:);

OutTableT.Time =  OutTableT.Time - OutTableT.Time(1);
time = OutTableT.Time;

tmp.Lin.V = 4:25;

for idx = 1 : length(tmp.Lin.V)
    timeIdx = time <= idx*100-10 & time >= idx*100-30;
    Lin_points.V(idx) = mean(OutTableT.Wind1VelX(timeIdx));

    Lin_points.RSpeed(idx) = mean(OutTableT.GenSpeed(timeIdx)/97/(60/(2*pi)));
    Lin_points.Pitch(idx) = mean(OutTableT.BlPitch1(timeIdx)/(180/(pi)));
    Lin_points.Torque(idx) = mean(OutTableT.GenTq(timeIdx));

end

% 
% Lin_points.RSpeed = max(Lin_points.RSpeed, Lin.RSpeed);
% Lin_points.Pitch = min(Lin_points.Pitch, Lin.Pitch);
% Lin_points.Torque = max(Lin_points.Torque, Lin.Torque);



figure;
subplot(4,1,1); plot(Lin_points.V,Lin_points.RSpeed,'bo', Lin.V, Lin.RSpeed,'*')
subplot(4,1,2); plot(Lin_points.V,Lin_points.Torque,'bo', Lin.V, Lin.Torque,'*')
subplot(4,1,3); plot(Lin_points.V,Lin_points.Pitch,'bo', Lin.V, Lin.Pitch,'*')
subplot(4,1,4); plot(Lin_points.V,Lin_points.Torque .* Lin_points.RSpeed ,'bo', Lin.V, Lin.Torque.*Lin.RSpeed ,'*')



if loadEqPoints ==  1
    Lin_points = Lin;
end

plotOn = 0;
Rotor_Lamda = ''; Rotor_Pitch = ''; Rotor_cQ = ''; Rotor_cT = '';
figDir = ''; titleOn = 0;

if useTurn == 0
    [wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, beta,Cq,Ct,Q3] = ....
        initModel5MWNREL(plotOn, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figDir,titleOn,multRb); %#ok<*ASGLU>

else
    [wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, pitch, Cq, Ct, Q3, Cp] = initModel5MWNRELTune(x);
end

K5 = [K,-K(:,4)];%#ok<*ASGLU>
DDNameCell = {'DD_test.sldd'; 'DD_Mdl1.sldd'};


for idxDD = 1: length(DDNameCell)
    DDName =  DDNameCell{idxDD}; %'DD_test.sldd'; %'DD_Mdl1.sldd'; %
    Q = Q3;
    DDDataNames = nameCell;

    myDictionaryObj = Simulink.data.dictionary.open(DDName);
    dDataSectObj = getSection(myDictionaryObj,'Design Data');

    for idx = 1: length(DDDataNames)
        eval(['tmp =',DDDataNames{idx},';']);
        try
            tempObj = getEntry(dDataSectObj,DDDataNames{idx});

            eval(['setValue(tempObj, ',DDDataNames{idx},')']);
        catch
            addEntry(dDataSectObj,DDDataNames{idx},tmp)
        end
    end

    saveChanges(myDictionaryObj)
end

% listEntry(myDictionaryObj)
%   Design Data   beta              DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   Ce                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   Cq                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   Ct                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   K5                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   lambda            DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   M                 DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   Q                 DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   rho               DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double
%   Design Data   wecs              DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          struct

% tmp = load('NREL5MW_linearised.mat');
% load('OutDataStep.mat');
% idxT = OutTable.Time >= 100;
% OutTableT = OutTable(idxT,:);
%
% OutTableT.Time =  OutTableT.Time - OutTableT.Time(1);
% time = OutTableT.Time;
%
% tmp.Lin.V = 4:25;
%
% for idx = 1 : length(tmp.Lin.V)
%     timeIdx = time <= idx*100-10 & time >= idx*100-30;
%     Lin_points1.V(idx) = mean(OutTableT.Wind1VelX(timeIdx));
%
%     Lin_points1.RSpeed(idx) = mean(OutTableT.GenSpeed(timeIdx)/97/(60/(2*pi)));
%     Lin_points1.Pitch(idx) = mean(OutTableT.BlPitch1(timeIdx)/(180/(pi)));
%     Lin_points1.Torque(idx) = mean(OutTableT.GenTq(timeIdx));
%
% end
%
% Lin.RSpeed = Lin_points1.RSpeed;
% Lin.Pitch = Lin_points1.Pitch;
% Lin.Torque = Lin_points1.Torque;
% Lin.V = Lin_points1.V;
%
% Lin_points = Lin;
