close all;
normStructNTW18 = runCompareModels(18,1);
close all;
normStructSweep = runCompareModels('Sweep',1);
close all;
runCompareCtrl('NTW18');
close all;
runCompareCtrl('Sweep');


%normStruct
% 
% normStruct = 
% 
%   struct with fields:
% 
%      RotSpeed: [0.0274 0.0317]
%        GenPwr: [0.0274 0.0318]
%         GenTq: [998.9999 998.9999]
%      BlPitch1: [0.1111 0.1271]
%     NcIMUTAxs: [2.0447 2.1530]
%     NcIMUTAys: [1.2814 1.2822]