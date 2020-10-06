clear; clc; close all

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

%% Dependen on index = 1
index = 1;
Vbar = Lin_points.V(index);
omegabar = Lin_points.RSpeed(index);
beta_bar = Lin_points.Pitch(index)/pi*180;
%Lambda_bar = omegabar*R/Vbar;