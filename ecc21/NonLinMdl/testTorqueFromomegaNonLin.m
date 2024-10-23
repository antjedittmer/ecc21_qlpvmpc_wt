%% Physical and turbine model constants
rho = 1.2250; % air density
N = 3;% Number of blades
Ng = 97; % Gearbox ratio
mh = 56780; % kg;  Hub mass
mb = 0.25*17740; %kg;  Modal mass of each blade
mn = 240000; % kg;  Nacellle mass
mtower = 347460;% kg; Mass of the tower and nacelle
mt = 0.25*mtower + mn + mh;
H =  87.6;

Jg = 534.116; %  kg*m^2; Inertia of the generator
Jr = 115926 + 3 * 11.776e6; % kg*m^2; Inertia of the rotor (Hub inertia + 3 blades)
f0 = 0.324; % Hz, First natural tower fore-aft frequency
% f0sw = 0.3120; % First natural tower sidewards frequency
wnb = 0.6993 * 2*pi; % rad/s First natural blade frequency 
wnt = f0 * 2*pi; % wnb; 0.3240
zetat = 1/100; %damping ratio of tower (Table 6.2)
zetab = 0.477465/100;% damping ratio of blade

Kt = wnt^2 * mt; % Stiffness of the tower s^2 + B/Ms + K/m
Bt = 2 * zetat * wnt *mt; % tower 2*6421;

Kb = wnb^2 * mb; %Stiffness of each blade
Bb = 2 * zetab* wnb *mb; %Damping of the blade
Ks = 867637000; %Nm/rad Stiffness of the transmission % 2*zeta*wn = B
Bs = 6215000;  %Nm/rad/sec %Damping of the transmission
rb = 63; % m blade radius
% etag = 0.944; %Drivetrain.Generator.Efficiency: 0.944;

% %% Effective wind speed and lambda
% Vbar = wind - zetadot*rb - ytdot;
% Lambda_bar = omegabar*rb/Vbar;
% 
% %% Grid points
% Cqb = interp2(lambdaVec,pitchVec,CqVec,Lambda_bar,pitch);
% Ctb = interp2(lambdaVec,pitchVec,CtVec,Lambda_bar,pitch);
% 
% %% 
% Tr = 0.5*rho*pi*rb^3*Cqb*Vbar^2;% br = dTr/domega
% Ft = 0.5*rho*pi*rb^2*Ctb*Vbar^2;% br = dTr/domega

%% Lagrange's Model matrices
% Force input w = [Ft_fa,Ft_sw,Tr,Tg]; Ft_fa = Ft, Ft_sw = 3/2*Tg
% States q: xdot_fa, zeta, xdot_sw, omega_r, omega_gr
M =[mt+N*mb N*mb*rb 0 0 0; % 
    N*mb*rb  N*mb*rb^2 0 0 0;
    0 0 mt 0 0;
    0 0 0 Jr 0;
    0 0 0 0 Jg*Ng^2];

Ce = [Bt 0 0 0 0;
    0 N*Bb*rb^2 0 0 0 ; %
    0 0 Bt 0 0
    0 0 0  Bs -Bs;
    0 0 0 -Bs Bs];

K = [Kt 0 0 0;
    0 N*Kb*rb^2 0 0; %
    0 0 Kt 0
    0 0 0 Ks;
    0 0 0 -Ks];

Q = [1 0 0 0;
    rb 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 -Ng];

L = [eye(4), [0;0;0;-1]];

A1 = [zeros(4), L;...
    -M\K, -M\Ce];
B1 = [zeros(4); M\Q];



