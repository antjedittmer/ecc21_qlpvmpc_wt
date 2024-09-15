function sys = createStateSpaceModel(Vbar,Rotor_Lamda,Rotor_Pitch,...
    Rotor_cQ,Rotor_cT,Cqdlambda,Ctdlambda,Cqdbeta,Ctdbeta, ...
    omegabar,Lambda_bar,beta_bar)
%createStateSpaceModel creates a 5DoF WT state space model at a given speed

%% Parameters for tower model
kfa = 1; %should this be a plus or a min?
Mfa = 4.38e5; % mt+N*mb: 400680
Cfa = 2*6421;
Kfa = 1.846e6;
kss = 0.1;
Mss = 4.38e5;
Css = 2*6421;
Kss = 1.846e6;

%% Define parameters
J = 4.2e7;% 3.5444E+07
R = 63;
rho = 1.225;
G = 97;
% Effic = 0.944;
Rotor_Lambda = Rotor_Lamda(1:150);
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
% Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');
%s = tf('s');

%% Derivatives
dlambdadomega = R/Vbar;                 %s
dlambdadV = -omegabar*R/Vbar^2;         %s/m
dlambdadxdotfa = -kfa*omegabar*R/Vbar^2; %s/m   xdotfa m/s

%% Grid points
dCqdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
dCtdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdlambda,Lambda_bar,beta_bar,'cubic');
dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
dCtdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic');
Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
Ctb = interp2(Rotor_Lambda,Rotor_Pitch,Ct,Lambda_bar,beta_bar,'cubic');

%% Calculation System 
kCT = 1/2*rho*pi*R^2;

A11=1/J*kCT*R*Vbar^2*dCqdlambdab*dlambdadomega;
A12=1/J*kCT*R*Vbar^2*dCqdlambdab*dlambdadxdotfa-2*kfa/J*kCT*R*Vbar*Cqb;

A21=1/Mfa*kCT*Vbar^2*dCtdlambdab*dlambdadomega;
A22=1/Mfa*kCT*Vbar^2*dCtdlambdab*dlambdadxdotfa-2*kfa/Mfa*kCT*Vbar*Ctb;

A31=kss/Mss*kCT*R*Vbar^2*dCqdlambdab*dlambdadomega;
A32=kss/Mss*kCT*R*Vbar^2*dCqdlambdab*dlambdadxdotfa-2*kss*kfa/Mss*kCT*R*Vbar*Cqb; % wrong

A=[A11 A12 0 0 0;
    A21 A22-Cfa/Mfa 0 -Kfa/Mfa 0;
    A31 A32  -Css/Mss 0 -Kss/Mss ;
    [0 1 0 0 0; 0 0 1 0 0]];

% Inputs: Torque [N/m], pitch angle beta [rad], wind speed V [m/s]
B11 = -G/J;
B12 = 1/J*kCT*R*Vbar^2*dCqdbetab;
B13 = 1/J*kCT*R*Vbar^2*dCqdlambdab*dlambdadV+2/J*kCT*R*Vbar*Cqb;
B22 = 1/Mfa*kCT*Vbar^2*dCtdbetab;
B23 = 1/Mfa*kCT*Vbar^2*dCtdlambdab*dlambdadV+2/Mfa*kCT*Vbar*Ctb;
B32 = kss/Mss*kCT*R*Vbar^2*dCqdbetab;
B33 = kss/Mss*kCT*R*Vbar^2*dCqdlambdab*dlambdadV+kss*2/Mss*kCT*R*Vbar*Cqb;

B=[B11 B12 B13;
    0   B22 B23;
    0   B32 B33;
    zeros(2,3);];

% Outputs: omega [rad/s], fore-aft dx_dot_f, fore-aft dx_dot_fa [m/s^2]
C(2,:)=A(2,:);
C(3,:)=A(3,:);
C(1,1)=1;

D(2,:)=B(2,:);
D(3,:)=B(3,:);

sys = ss(A,B,C,D);