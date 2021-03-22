function [currentOutTableTest, tictoc_LPVMPC ]= getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname1,varnames)
% OutTableTest1 retrieves simulation data, either from previous runs or by
% running the simulation. All inputs are required.
% - matFileOutTableTest1: mat file name
% - loadData: data loaded if available
% - OutTable: FAST simulation table
% - simMdlname1: Simulink mdl name
% - varnames: variables names of simulation output

% Get the name of the table
[~, strTable] = fileparts(matFileOutTableTest1);

% Get the simulation output
tictoc_LPVMPC = [];
if exist(matFileOutTableTest1,'file') == 2 && loadData
    tmp = load(matFileOutTableTest1,strTable);
    currentOutTableTest = tmp.(strTable);
    if isfield(tmp.tictoc_LPVMPC) 
        tictoc_LPVMPC = tmp.tictoc_LPVMPC;
    end
else
    WS = 'base';
    % assignin(WS,'DT',DT)
    assignin(WS,'OutTable',OutTable)
    open_system(simMdlname1);
    out1 = sim(simMdlname1,'SrcWorkspace',WS);
    currentOutTableTest = array2table(out1.OutDataTest,'VariableNames',varnames);
    eval([strTable,' = currentOutTableTest;']);
    if contains(strTable,'MPC')
        save(matFileOutTableTest1,strTable,'tictoc_LPVMPC');
    else
    save(matFileOutTableTest1,strTable);
    end
    
end
