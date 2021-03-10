function [wecs,rho, tau, kappa ] = initModel5MWNREL

% Set path to inputdata directory
mainDir = fullfile(pwd);
dataInDir = fullfile(fileparts(pwd),'dataIn');
addpath(dataInDir);

% Set path to figure directory
%figDir = fileparts(mainDir); %'html','figDir');
figDir = fullfile(fileparts(mainDir),'figDir');

if ~isfolder(figDir)
    mkdir(figDir)
end

load('NREL5MW_CPdata_new3','Lin_points','Rotor_Lamda','Rotor_Pitch','Rotor_cP','Rotor_cQ','Rotor_cT');
%load('NREL5MW_CPdata','Rotor_Lamda','Rotor_Pitch','Rotor_cP','Rotor_cQ','Rotor_cT');

%% Define parameters
% 5MW FAST WIND TURBINE (NREL/TP-500-38060) Definition of a 5-MW Reference Wind Turbine for Offshore
% System Development
% Cut-In, Rated, Cut-Out Wind Speed 3 m/s, 11.4 m/s, 25 m/s
% Cut-In, Rated Rotor Speed 6.9 rpm, 12.1 rpm
%wr,wg in rad/s, Tg in kN, Pg in kW,
wecs.Ng = 97; %Gearbox ratio
wecs.mh =  56780; % hub mass
wecs.mb =  0.25*17740; %kg Modal mass of each blade
wecs.mn = 240000; % Nacellle mass
wecs.mtower = 347460;%347460;%kg Mass of the tower and nacelle
wecs.mt = 0.25*wecs.mtower+wecs.mn+wecs.mh;
wecs.H =  87.6;

wecs.Jrg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator
wecs.Jr = 115926 + 3 * 11.776e6 + wecs.Jrg * wecs.Ng^2; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
f0 = 0.324; %First natural tower frequency
% f0sw = 0.3120;
wecs.wnb = 1 * 2*pi; %2.0205; %Natural frequency of blade 0.6993
wecs.wnt = f0 * 2*pi; % wecs.wnb; 0.3240
wecs.zetat = 1/100; %damping ratio of tower (Table 6.2)
wecs.zetab = 0.477465/100;% damping ratio of blade

wecs.Kt = wecs.wnt^2 * wecs.mt; % Stiffness of the tower s^2 + B/Ms + K/m
wecs.Bt = 2 *wecs.zetat * wecs.wnt *wecs.mt; % tower 2*6421;
wecs.Kb = wecs.wnb^2 * wecs.mb; %Stiffness of each blade
wecs.Bb = 2 *wecs.zetab*wecs.wnb *wecs.mb ; %Damping of the blade
wecs.Ks = 867637000; %Nm/rad Stiffness of the transmission
% 2*zeta*wn=B
wecs.Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
wecs.N = 3;%Number of blades
wecs.rb = 63; % m blade radius

rho = 1.225; % Air density (kg/m^3)
% Bg = 0.9; % Tg= Bg(wg - wz). Unused because we use Tg as input
tau = 0.1; % for pitch actuator
kappa = 0.01;

% Tower equivalent mass, MT 438,000 kg
% Tower equivalent damping, CT 6421 wecs.N s/m
% Tower equivalent stiffness, KT 1,846,000 wecs.N/m

%% Lagrange's Model matrices
% Force input w = [Ft,3/2*Tg, Tr,Tg];
% States q: xdot_fa, zeta, xdot_sw, omega_r, omega_gr
% States
M =[wecs.mt+wecs.N*wecs.mb wecs.N*wecs.mb*wecs.rb 0 0 0; % M(1,1) mt+wecs.N*mb  -> mt (already includes blades)
    wecs.N*wecs.mb*wecs.rb  wecs.N*wecs.mb*wecs.rb^2 0 0 0;
    0 0 wecs.mt 0 0;
    0 0 0 wecs.Jr 0;
    0 0 0 0 wecs.Jrg*wecs.Ng^2];

Ce = [wecs.Bt 0 0 0 0;
    0 wecs.N*wecs.Bb*wecs.rb^2 0 0 0 ; %
    0 0 wecs.Bt 0 0
    0 0 0  wecs.Bs -wecs.Bs;
    0 0 0 -wecs.Bs wecs.Bs];

K = [wecs.Kt 0 0 0;
    0 wecs.N*wecs.Kb*wecs.rb^2 0 0; %
    0 0 wecs.Kt 0
    0 0 0 wecs.Ks;
    0 0 0 -wecs.Ks];

Q = [1 0 0 0;
    wecs.rb 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 -wecs.Ng];

L = [eye(3), [0;0;-1]];

A1 = [zeros(3), L;...
    -M\K, -M\Ce];
B1 = [zeros(3); M\Q];


tau = 0.1; % for pitch actuator
tau_g  = 0.01;

Drivetrain.Generator.Efficiency = 0.944;


%% Model
Rotor_Lambda = Rotor_Lamda(1:150);
idxPitch = Rotor_Pitch>= -8;
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch(idxPitch));
Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');

% % Bindu: Calculate CP matrix for interpolated surf plot
% Cp2Comp = interp2(Rotor_lambda2,Rotor_beta2,Rotor_cP2,Xq,Yq,'nearest');
% Cq2Comp = interp2(Rotor_lambda2,Rotor_beta2,Rotor_cQ2,Xq,Yq,'nearest');
% Ct2Comp = interp2(Rotor_lambda2,Rotor_beta2,Rotor_cT2,Xq,Yq,'nearest');

% xlabelStr = 'tip speed ratio \lambda [-]';
% ylabelStr = 'pitch angle \beta [grad]';
% zlabelStr = 'Force coefficient C_T [-]';
% titleStrData1Interp = ['Aerodynamique ',zlabelStr,'  from NREL FAST 5 MW'];
% 
% figure(1); surf(Xq,Yq,Ct,'FaceColor','interp','EdgeColor','none'); hold on
% axis tight;
% xlabel(xlabelStr); ylabel(ylabelStr ); zlabel(zlabelStr);
% title(titleStrData1Interp);
% view(60,30)
% 
% print(gcf,[fullfile(figDir,'Ct_NRRLFAST5MW')], '-dpng');
% print(gcf,[fullfile(figDir,'Ct_NRRLFAST5MW')], '-depsc');
% 
% 
% zlabelStr = 'Torque coefficient C_Q [-]';
% titleStrData1Interp = ['Aerodynamique ',zlabelStr,' from NREL FAST 5 MW'];
% 
% figure(2); surf(Xq,Yq,Cq,'FaceColor','interp','EdgeColor','none'); hold on
% axis tight;
% xlabel(xlabelStr); ylabel(ylabelStr ); zlabel(zlabelStr);
% title(titleStrData1Interp)
% view(60,30)
% 
% print(gcf,[fullfile(figDir,'Cq_NRRLFAST5MW')], '-dpng');
% print(gcf,[fullfile(figDir,'Cq_NRRLFAST5MW')], '-depsc');


%% Compute gradients
[Cqdlambda,Cqdbeta] = gradient(Cq,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Ctdlambda,Ctdbeta] = gradient(Ct,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
Clambda = Cq./Xq.^2;
[Clambdadlambda,Clambdadpitch] = gradient(Clambda,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
Rotor_Lamda = Rotor_Lambda;


%% Some calculated parameters from FAST
% 
%     Hub-Height            (m)            90.437
%     Flexible Tower Length (m)            87.600
%     Flexible Blade Length (m)            61.500
% 
% 
% Rotor mass properties:
% 
%     Rotor Mass            (kg)       109992.734
%     Rotor Inertia         (kg-m^2) 39118660.000
%                                         Blade 1      Blade 2      Blade 3
%                                         -------      -------      -------
%     Mass                  (kg)        17737.578    17737.578    17737.578
%     Second Mass Moment    (kg-m^2) 11886951.000 11886951.000 11886951.000
%     First Mass Moment     (kg-m)     366280.750   366280.750   366280.750
%     Center of Mass        (m)            20.650       20.650       20.650
% 
% 
% Additional mass properties:
% 
%     Tower-top Mass        (kg)       349992.750
%     Tower Mass            (kg)       348398.125
%     Platform Mass         (kg)            0.000
%     Mass Incl. Platform   (kg)       698390.875

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


%% Unused code: Parameters for 2nd order tower model from TU Delft
% kfa = 1; 
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

