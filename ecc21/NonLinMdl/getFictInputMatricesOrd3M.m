function [B1,A1,uLin,Aout,Bout] = getFictInputMatricesOrd3M(wind,Tg,pitch,lambda,beta,Cq,Ct,DT,omegabar,ytdot,zetadot)

if ~nargin
    LinPointsMFile;
MdlMFile;
idx = 15;
wind = Lin_points.V(idx);
Tg = Lin_points.Torque(idx)*10^3;
pitch = Lin_points.Pitch(idx);
DT = 0.008;
omegabar = Lin_points.RSpeed(idx);
ytdot = 0;
zetadot = 0;

end


[B1,A1,uLin,Aout,Bout] = getFictSub(wind,Tg,pitch,lambda,beta,Cq,Ct,DT,omegabar,ytdot,zetadot);

function [B1,A1,uLin,Aout,Bout] = getFictSub(wind,Tg,pitch,lambda,beta,Cq,Ct,DT,omegabar,ytdot,zetadot)
% getTrFt
% representation for ve
%
% Inputs: * Ts: Sampling time
%         * AwithV: Scheduling parameter vector: AwithV = V
%
% Outputs: * A: System matrix A(AwithV) at current estimated air speed AwithV, discretised with Ts
%          * B: Input matrix B(AwithV) at current AwithV, discretised with Ts

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for WECS simulation: Antje Dittmer


persistent lambdaVec pitchVec CqVec Cqdlambda Cqdbeta CtVec

if isempty(lambdaVec)
    lambdaVec =  lambda;
    pitchVec = beta;
    CqVec = Cq;
    CtVec = Ct;
    %% Compute gradients %%[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
    [Cqdlambda,Cqdbeta] = gradient(CqVec,mean(diff(lambda)),mean(diff(beta))); %dcQ/dbetaDeg * dbetaDeg/dbetaRad *betaRAd

end

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
Jgr = Jg * Ng^2; % Inertia of the generator in LSS
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

%% Effective wind speed and lambda
Vbar = wind - zetadot*rb - ytdot;
Lambda_bar = omegabar*rb/Vbar;
dlambdadomega = rb/Vbar;
dlambdadV = - omegabar*rb/Vbar^2;
beta_bar = pitch;

%% Grid points
Cqb = interp2(lambdaVec,pitchVec,CqVec,Lambda_bar,beta_bar);
dCqdlambdab = interp2(lambdaVec,pitchVec,Cqdlambda,Lambda_bar,beta_bar);
dCqdbetab = interp2(lambdaVec,pitchVec,Cqdbeta,Lambda_bar,beta_bar);
Tr = 0.5*rho*pi*rb^3*Cqb*Vbar^2;% br = dTr/domega
%Ft = 0.5*rho*pi*rb^2*Ctb*Vbar^2;% br = dTr/domega

%% Tr=0.5*pi*rb^3*Cq*v^2
dTrdomega = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;% br = dTr/domega
dTrdbeta = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;% Krb = dTr/dbeta
dTrdV =  2*0.5*rho*pi*rb^3*Vbar*Cqb + ...
    0.5*rho*pi*rb^3*Vbar^2*dCqdlambdab*dlambdadV;

%% Lagrange's Model matrices
% Force input w = [Ft_fa,Ft_sw,Tr,Tg]; Ft_fa = Ft, Ft_sw = 3/2*Tg
% States q: xdot_fa, zeta, xdot_sw, omega_r, omega_gr
vecOmega = 4:5; 
Mall =[mt+N*mb N*mb*rb 0 0 0; % 
    N*mb*rb  N*mb*rb^2 0 0 0;
    0 0 mt 0 0;
    0 0 0 Jr 0;
    0 0 0 0 Jg*Ng^2];
M = Mall(vecOmega,vecOmega);


Ceall = [Bt 0 0 0 0;
    0 N*Bb*rb^2 0 0 0 ; %
    0 0 Bt 0 0
    0 0 0  Bs -Bs;
    0 0 0 -Bs Bs];
Ce = Ceall(vecOmega,vecOmega);

Kall = [Kt 0 0 0;
    0 N*Kb*rb^2 0 0; %
    0 0 Kt 0
    0 0 0 Ks;
    0 0 0 -Ks];
K = Kall(vecOmega,vecOmega(1));


Qall = [1 0 0 0;
    rb 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 -Ng];
Q = Qall(vecOmega,vecOmega-1);

Lall = [eye(4), [0;0;0;-1]];
L = Lall(vecOmega(1),vecOmega);
%L(vecOmega,vecOmega(1)) = Lall;

nqt = size(L,1);
A1 = [zeros(nqt), L;...
    -M\K, -M\Ce];
B1 = [zeros(size(L)); M\Q];

%% Calculate discrete forms of A,B
% A2 = A1*A1;
% A3 = A1*A2;
% Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*Arho+DT*eye(size(Arho));
% %% Finish calculation of A,B
% Ak = Apre*A1 + eye(size(A1));
% Bk = Apre*B1;
uLin = [Tr,Tg]';

%% Create system matrix A:
% States: [Theta_s [rad], omega_r [rad/s], omega_gr [rad/s], T_g [Nm]]
% States for model xA = [y, xdot] with y = [Tg, beta]
% x = [ThetaS, OmegaR, OmegaG, Tg, beta,]
Atemp = [[0,1,-1];
    [-Ks/Jr, (-Bs + dTrdomega)/Jr, Bs/Jr,];
    Ks/Jgr, Bs/Jgr, -Bs/Jgr];

%% Create input matrix B
% Input orig model: V [m/s], T_g_ref [Nm], beta [rad],
% -> vel model: Vdot [m/s^2]  , T_g_refDot [Nm/s]
% for beta: Vdot [m/s^2]  T_g_refDot [Nm/s] betaDot [rad/s],
Btemp = [zeros(1,3); dTrdV/Jr, 0 dTrdbeta/Jr;  0,  -Ng/Jgr,0];

%% Check rank
Bu = Btemp(:,2:3);
Bd = Btemp(:,1);
nxEval = size(Atemp,2);

% Also create an additional state for the disturbance
tau1 =  1/DT;
Aout = [[Atemp, zeros(nxEval,1)]; [zeros(1,nxEval),-tau1]];
E = [Bd;tau1]; %xd_dot = -xd + Ed
B = [Bu; zeros(1,2)];
Bout = [E,B];