% clear; clc; close all
% addpath('C:\Program Files\export_fig'); %C:\Users\bs\Desktop\export_fig-master\export_fig-master;
% addpath C:\Werk_Data\Software\hline_vline;
% load 'NREL5MW_linearised_reduced'
load('NREL5MW_CPdata_new3','Rotor_Lamda','Rotor_Pitch','Rotor_cP','Rotor_cQ','Rotor_cT','Lin_points');
load('OutDataSweep.mat');

%% Define parameters
% 5MW FAST WIND TURBINE (NREL/TP-500-38060)
% Definition of a 5-MW Reference Wind Turbine for Offshore System Development
% Cut-In, Rated, Cut-Out Wind Speed 3 m/s, 11.4 m/s, 25 m/s
% Cut-In, Rated Rotor Speed 6.9 rpm, 12.1 rpm
%wr,wg in rad/s, Tg in kN, Pg in kW,
Ng = 97; %Gearbox ratio
mt = 347460 +240000;%347460;%kg Mass of the tower and nacelle
mb = 17740; %kg Mass of each blade
Jr = 115926 + 3 * 11.776e6; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator
wnt = 2.1; %2.1 ; %2.9003; %Natural frequency of tower
wnb = 2.0205; %Natural frequency of blade
zetat = 0.01; %damping ratio of tower
zetab = 0.00477465;%damping ratio of blade
% K=wn^2
Kt = wnt^2 * mt; % Stiffness of the tower s^2 + B/Ms + K/m
Kb = wnb^2 * mb; %Stiffness of each blade

Ks = 867637000; %Nm/rad Stiffness of the transmission
% 2*zeta*wn=B
Bt = 2*zetat*wnt *mt; % Damping of the tower
Bb = 2*zetab*wnb *mb ; %Damping of the blade
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission

N = 3;%Number of blades
rb = 63; % m blade radius
R = rb;
rho = 1.225; % Air density (kg/m^3)
Bg = 0.9; % Tg= Bg(wg - wz).
tau = 0.1; % for pitch actuator
tau_g  = 0.01;

Drivetrain.Generator.Efficiency = 0.944;


%% Model
M= [mt+N*mb N*mb*rb 0 0;
    N*mb*rb N*mb*rb^2 0 0;
    0 0 Jr 0;
    0 0 0 Jg*Ng^2]; % 

Ce = [Bt 0 0 0;
    0 N*Bb*rb^2 0 0;
    0 0 Bs -Bs;
    0 0 -Bs Bs];

K = [Kt 0 0 0;
    0 N*Kb*rb^2 0 0;
    0 0 Ks -Ks;
    0 0 -Ks Ks];

Q = [1 0 0;    
    R 0 0;
    0 1 0
    0 0 -Ng];


tildaK = [Kt 0 0;
    0 Kb*rb^2 0;
    0 0 Ks;
    0 0 -Ks];

M34 = M(3:4,3:4);
K34 = K(3:4,3:4);
Ce34 = Ce(3:4,3:4);
Q34 = Q(3:4,2:3);

Rotor_Lambda = Rotor_Lamda(1:150);
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');


%% Compute gradients
[Cqdlambda,Cqdbeta] = gradient(Cq,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Ctdlambda,Ctdbeta] = gradient(Ct,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
Clambda = Cq./Xq.^2;
[Clambdadlambda,Clambdadpitch] = gradient(Clambda,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
Rotor_Lamda = Rotor_Lambda;

OutTable = array2table(OutData,'VariableNames', OutList);

idxTime = 1: height(OutTable);
idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

%% Unused code: For tests
% CEM = M\Ce;
% KM = M\K;
% QM = M\Q;
% 
% Atest = [[zeros(2), eye(2)]; [-KM(3:4,3:4), -CEM(3:4,3:4)]];
% Btest = QM(:,2);
% 
% Ctest = zeros(1,4);
% Ctest(3) = 1;
% Dtest = zeros(1,1);


%% Unused code:Parameters for tower model
% kfa = 1; %should this be a plus or a min?
% Mfa = 4.38e5; % mt+N*mb: 400680
% Cfa = 2*6421;
% Kfa = 1.846e6;
% 
% % mb =(Mfa - mt)/N;
% mt = Mfa;
% Kt = Kfa;
% Bt = Cfa; % 2*zetat*wnt *mt

% Kb = wnb^2*mb; %
% Bb  = 2*zetab*wnb*mb;

