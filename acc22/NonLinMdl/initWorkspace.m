%% Clear workspace
clear; close all; clc; 
restoredefaultpath;
bdclose all;
clear ABgL_WT; % LUT variables saved as persistent variables
Simulink.data.dictionary.closeAll;

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
dirSub = {'subfunctionsMPC', 'subfunctionsUtil', 'WECSCtrl', 'WECSMDl','WECSAct'}; 

% Get relative paths with respect to testCrl.m file
dirFileName = mfilename('fullpath');
dirWork = fileparts(dirFileName);
dirParent = fileparts(dirWork);

% Add parentDir and workdir to path: Can be called in both directories
addpath(dirParent)
addpath(dirWork)

% Add LinMdl (for WECS param in workspace for testing)
addpath(fullfile(dirParent,'LinMdl'));

% Add path to input data (extracted from FAST with FASTtool)
dataInPath = fullfile(dirParent,dataIn);
addpath(dataInPath);

% Add output path to save mat files
dataOutPath = fullfile(dirWork,dataOut); 
if ~isdir(dataOutPath) %#ok<ISDIR>
    mkdir(dataOutPath);
end
addpath(dataOutPath);

% Add subdirerectoreis
for dirIdx = 1: length(dirSub)
    addpath(fullfile(dirWork,dirSub{dirIdx}));
end

%% Load FASTTool wind as test signal
DT = 0.008; % one sample time 
% load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');
load(fullfile(dataInPath,'OutDataWind18NTWnoise.mat'),'OutTable');

%% Load variables in workspace (this might go to a data dictionary)
LinPointsMFile;
MdlMFile;
ControllerMFile;
switchNoise = 1;

%% Create UIO
Jgr = wecs.Jg * wecs.Ng^2;

A = [0,1,-1 ;
    -wecs.Ks/wecs.Jr, -wecs.Bs/wecs.Jr , wecs.Bs/wecs.Jr;
    wecs.Ks/Jgr, wecs.Bs/Jgr, -wecs.Bs/Jgr];
B = [0; 0; -wecs.Ng/Jgr];
E = [0; 1/wecs.Jr;0 ];
C = [0 1 0;0 0 1]; % 
%fprintf('rank(C*E) is %d; rank(E) is %d\n', rank(C*E),rank(E));

H = E*inv((C*E)'*(C*E))*(C*E)'; %#ok<MINV,*NOPTS> display output, inv is fast enough
nx = size(A,1);
T = eye(nx)-H*C;
A1 = T*A;

% Pole placement with Q and R
Qe = eye(nx); %diag([1,10^11,10^9])*[B,E]*[B,E]'; %eye(nx); %
% r1 = 10^-2;
% r2 = 10^-2;
Re = eye(2); %diag([r1,r2]); %*100;
K1 = lqr(A',C',Qe,Re)';


% Calculate F and K
F = A1 - K1*C; %T*A -K1*C
K = K1 + F*H; %T = eye(size(H*C)) - H*C;

% Sys with no ctrl input: xdot = Ax + Ed, y = Cx
nd = size(E,2); %nu = 0;
nu = size(B,2);
ny = size(C,1);
sys = ss(A,[E,B],C,zeros(ny,nd+nu)); %c2d(ss(A,[E,B],C,zeros(ny,nd+nu)),DT,'zoh'); % ss(A,[E,B],C,zeros(ny,nd+2));

% Est with with control input: zdot = Fz + TBu + Ky, xhat = z + Hy
sysEstTr = c2d(ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]),DT,'zoh'); %ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]); %


%% Clean-up
% clear RESTORE*
% clear data*
% clear dir*



