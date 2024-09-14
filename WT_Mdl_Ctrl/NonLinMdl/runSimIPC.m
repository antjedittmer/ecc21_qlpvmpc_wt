function [Won,Kon,refOn,OutTableON] = runSimIPC(xCtrl,ref,OutList,weightVec)

if nargin < 4
    weightVec = [1.25,1,0.5,1,1.25];
end

%IPC_K_1P = xCtrl(1); %#ok<*NASGU>
IPC_KI_1p = xCtrl(1);%#ok<*NASGU>
IPC_KP_1p = xCtrl(2);
kp_Turm = xCtrl(3);
% IPC_KI_1p_tilt = xCtrl(5);
% IPC_KP_1p_tilt = xCtrl(6);

hEntry = who;
for idx = 1: length(hEntry)
    assignin('base',hEntry{idx},eval(hEntry{idx}))
end

WS = 'base';
out1 = sim('test_SimulinkMdl3LPV_IPCa','SrcWorkspace',WS);
 
OutTableON = array2table(out1.OutDataLPV,'VariableNames',OutList);

[Won,Kon,refOn] = getEvalCriteria(OutTableON,ref,weightVec);
