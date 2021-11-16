
%https://onlinelibrary.wiley.com/doi/epdf/10.1002/we.1533
Ng = 97; %Gearbox ratio
mt = 347460 +240000;%347460;%kg Mass of the tower and nacelle
mb = 17740; %kg Mass of each blade
Jr = 115926 + 3 * 11776047;%3 * 11.776e6; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator
wnt = 2*pi*.32;%2.1; %2.1 ; %2.9003; %Natural frequency of tower
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


%% TU Delf
%J=115926 + 3 * 11776047+534.116*97*97=4.047e+07 accprding to paper
J=4.2e7;