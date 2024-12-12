% script to update data dictionary
% clear;
function updateDDMdl1(multRb)

if ~nargin
    multRb = 0.75;
end

plotOn = 0;
Rotor_Lamda = ''; Rotor_Pitch = ''; Rotor_cQ = ''; Rotor_cT = '';
figDir = ''; titleOn = 0;

[wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, beta,Cq,Ct,Q3] = ....
    initModel5MWNREL(plotOn, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figDir,titleOn,multRb); %#ok<*ASGLU> 
K5 = [K,-K(:,4)];%#ok<*ASGLU> 
DDNameCell = {'DD_test.sldd'; 'DD_Mdl1.sldd'};

% tmp = load('NREL5MW_linearised.mat');
% Lin_points = tmp.Lin;

for idxDD = 1: length(DDNameCell)
    DDName =  DDNameCell{idxDD}; %'DD_test.sldd'; %'DD_Mdl1.sldd'; %
    Q = Q3;
    DDDataNames = {'beta';'Ce';'Cq';'Ct';'K5'; 'lambda'; 'M';'Q'; 'rho';'wecs'; 'Lin_points'};



    myDictionaryObj = Simulink.data.dictionary.open(DDName);
    dDataSectObj = getSection(myDictionaryObj,'Design Data');

    for idx = 1: length(DDDataNames)
        try
            tempObj = getEntry(dDataSectObj,DDDataNames{idx});

            eval(['setValue(tempObj, ',DDDataNames{idx},')']);
        catch
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


