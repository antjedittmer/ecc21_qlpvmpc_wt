% 'omega_r, omega_g, theta_dt'
windspeed = 23;
R = 33.29; %meter
rho = 1.225; %Air density
area = pi*R^2;%3481.588; % Area 

%% drive train
Kdt = 31.8*10^6; %Nm/rad
Dampdt = 212.2;
Jr = 1.86*10^6; % kgm^2
Jg = 56.29; % kgm^2
Ng = 74.38;

energyinitial = 0.5*Jg*0.01^2 + 0.5*Jr*0.01^2 ;

Adt = [0 0 -Kdt/Jr;
        0 0 Kdt/(Ng*Jg);
        1 -1/Ng 0];
Bdt = [1/Jr 0;
        0 -1/Jg;
        0 0];    
Cdt = eye(3);
Ddt = zeros(3,2);
dtinitial = [3.04;(3.04*74.38);0];
%dtinitial = [0.01;0.01;0];

% Adt = [-Dampdt/Jr Dampdt/(Ng*Jr) -Kdt/Jr;
%         Dampdt/(Ng*Jg) -Dampdt/(Ng^2*Jg) Kdt/Jg;
%         1 -1/Ng 0];
% Bdt = [1/Jr 0;
%         0 -1/Jg;
%         0 0];    
% Cdt = eye(3);
% Ddt = zeros(3,2);
% dtinitial = [0.1;0.1;0];

% Adt = 0*eye(2);
% Bdt = [1/Jr 0;
%         0 -1/Jg];
% Cdt = eye(2);
% Ddt = 0*eye(2);
% dtinitial = [0.1;0.1];
%%
eta = 1;
tau_g = 0.02; %sec
delT=0.05;

%% tower dynamics
Kt=1652000;%3.5159e3;
Mt = 54606;%4.2278e5;
Bt=7589;%2.0213e3;
hubheight = 50; %m

Atw = [0 1;
    -Kt/Mt -Bt/Mt];
Btw = [0; 1/Mt];
Ctw = eye(2);
Dtw = zeros(2,1);

x_0_tw = [0.05 0]';

%% faults
biasfault = 10;
Constant_tau_gc = 2000;