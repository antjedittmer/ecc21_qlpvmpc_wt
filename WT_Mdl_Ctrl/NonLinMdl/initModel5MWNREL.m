% clear; clc; close all

% load 'NREL5MW_linearised_reduced'
% load 'NREL5MW_linearised'
load('NREL5MW_CPdata_new3','Rotor_Lamda','Rotor_Pitch','Rotor_cP','Rotor_cQ','Rotor_cT','Lin_points');
% Lin_points = Lin;

%% Define parameters
J = 4.2e7;% 3.5444E+07
R = 63;
rho = 1.225;
G = 97;
Effic = 0.944;
H = 87.6; %height above ground
Rotor_Lambda = Rotor_Lamda(1:150);
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');
%s = tf('s');

%% Compute gradients
[Cqdlambda,Cqdbeta] = gradient(Cq,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Ctdlambda,Ctdbeta] = gradient(Ct,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
Clambda=Cq./Xq.^2;
[Clambdadlambda,Clambdadpitch] = gradient(Clambda,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));


%% Parameters for tower model
kfa = 1; %should this be a plus or a min?
Mfa = 4.38e5; % mt+N*mb: 400680
Cfa = 2*6421;
Kfa = 1.846e6;
kss = 0.1;
Mss = 4.38e5;
Css = 2*6421;
Kss = 1.846e6;

% ---------------------- MASS AND INERTIA ----------------------------------------
%           0   TipMass(1)  - Tip-brake mass, blade 1 (kg)
%           0   TipMass(2)  - Tip-brake mass, blade 2 (kg)
%           0   TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]
%  5.678E+04      HubMass     - Hub mass (kg)
%  1.159E+05      HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)
%  534.116      GenIner     - Generator inertia about HSS (kg m^2)
%  2.400E+05      NacMass     - Nacelle mass (kg)
%  2.599E+06      NacYIner    - Nacelle inertia about yaw axis (kg m^2)
%           0   YawBrMass   - Yaw bearing mass (kg)
%           0   PtfmMass    - Platform mass (kg)
%           0   PtfmRIner   - Platform inertia for roll tilt rotation about the platform CM (kg m^2)
%           0   PtfmPIner   - Platform inertia for pitch tilt rotation about the platform CM (kg m^2)
%           0   PtfmYIner   - Platform inertia for yaw rotation about the platform CM (kg m^2)
% 

%% Dependen on index = 1
index = 1;
Vbar = Lin_points.V(index);
omegabar = Lin_points.RSpeed(index);
beta_bar = Lin_points.Pitch(index)/pi*180;
%Lambda_bar = omegabar*R/Vbar;