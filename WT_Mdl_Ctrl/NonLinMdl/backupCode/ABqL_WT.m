function [A,B] = ABqL_WT(DT,rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ)
% ABqL calculates system (A) and input (B) matrix for quasiLPV
% representation for ve
%
% Inputs: * Ts: Sampling time
%         * AwithV: Scheduling parameter vector: AwithV = V
%
% Outputs: * A: System matrix A(AwithV) at current estimated air speed AwithV, discretised with Ts
%          * B: Input matrix B(AwithV) at current AwithV, discretised with Ts

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%% Parameters
Ng = 97; %Gearbox ratio

if length(rho) == 1
    Vbar = rho;
    
    % LUT table rows
    LUT_V = 4:25;
    LUT_omega = [0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)];
    LUT_beta = [zeros(1,8),0.0642    0.1126    0.1493    0.1808    0.2093    0.2356    0.2602...
        0.2833    0.3054    0.3269    0.3477    0.3676    0.3868    0.4056];  
    
    omegabar = interp1(LUT_V,LUT_omega,Vbar,'pchip');
    beta_bar = interp1(LUT_V,LUT_beta,Vbar,'pchip')/pi*180;
    
elseif length(rho) == 3
    % rho = xdotk_1([1,3,4]); % Tg not included: [Omega [1/s], beta[rad], V [m/s]];
    Vbar = rho(3);
    omegabar = rho(1)/Ng; % OmegaG/Ng as approximation of omegaR
    beta_bar = rho(2)/pi*180;
    
else
    error('Unexpected length of scheduling parameter rho')
end

  
rb = 63; % m blade radius
Lambda_bar = omegabar*rb/Vbar;
Ng = 97; %Gearbox ratio
% mt = 347460 +240000;%347460;%kg Mass of the tower and nacelle
% mb = 17740; %kg Mass of each blade
Jr = 115926 + 3 * 11.776e6; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
Jgr = Jg * Ng^2; % Inertia of the generator in LSS
% wnt = 2.1; %2.1 ; %2.9003; %Natural frequency of tower
% wnb = 2.0205; %Natural frequency of blade
% zetat = 0.01; %damping ratio of tower
% zetab = 0.00477465;%damping ratio of blade
% Kt = wnt^2 * mt; % Stiffness of tower s^2 + B/Ms + K/m; Needed for tower
% Kb = wnb^2 * mb; %Stiffness of each blade
Ks = 867637000; %Nm/rad Stiffness of the transmission
% 2*zeta*wn=B
% Bt = 2*zetat*wnt *mt; % Damping of the tower
% Bb = 2*zetab*wnb *mb ; %Damping of the blade
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission

% N = 3;%Number of blades
rb = 63; % m blade radius
rho = 1.225; % Air density (kg/m^3)
%Bg = 0.9; % Tg= Bg(wg - wz).
tau_g = 0.01; % for pitch actuator and airspeed
tau = 0.01;

%%
dlambdadomega= rb/Vbar;
dlambdadV= - omegabar*rb/Vbar^2;

[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
% Cp = interp2(Rotor_Lambda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lambda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
% Ct = interp2(Rotor_Lambda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');


%% Compute gradients %%[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Cqdlambda,Cqdbeta] = gradient(Cq,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
% [Ctdlambda,Ctdbeta] = gradient(Ct,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
% Clambda = Cq./Xq.^2;
% [Clambdadlambda,Clambdadpitch] = gradient(Clambda,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));


%% Grid points
dCqdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
%dCtdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic');
Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');


%% Tr=0.5*pi*rb^3*Cq*v^2
dTrdomega = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;% br = dTr/domega
dTrdV = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadV*Vbar^2+2*Vbar*Cqb);% Krv =  dTr/dV
dTrdbeta = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;%Krb = dTr/dbeta


%% Create system matrix A:
% States: [Theta_s [rad], omega_r [rad/s], omega_g [rad/s], T_g [Nm]]
% States for model
% xA = [y, xdot] with
% y = [omegaG, Tg, beta, V]
% x = [ThetaS, OmegaR, OmegaG, Tg, beta, V]
% xdot = [ThetaS, OmegaR, OmegaG, Tg, beta, V]';
% u = [V,beta_ref,TgRef]; Input: V [m/s], beta [rad], T_g_ref [Nm] (in HSS)

Aeval = [[0,1,-1,0] ;
    [-Ks/Jr, (-Bs + dTrdomega)/Jr, Bs/Jr,0];
    Ks/Jgr, Bs/Jgr, -Bs/Jgr, -Ng/Jgr;
    zeros(1,3), -1/tau_g];


%% Create input matrix B
% Input orig model: V [m/s], beta [rad], T_g_ref [Nm]
% Input vel model: Vdot [m/s^2]  betaDot [rad/s], T_g_refDot [Nm/s]
% Assumption Vdot = 0
% Output: omega_g [rad/s], T_g [Nm]
% For orig model:
% Beval = [zeros(1,3);
%     dTrdV/Jr, dTrdbeta/Jr,0;
%     zeros(1,3);
%     0, 0, 1/tau_g];

Beval = [zeros(1,2);
    dTrdbeta/Jr,0;
    zeros(1,2);
    0, 1/tau_g];


%% Add airspeed V as a state with Vdot = -1/tau V+ 1/tau V
% Add beta as state
% AwithBeta = [Aeval, Beval(:,1); [zeros(1,size(Aeval,2)),-1/tau]]; %
%       
% %Inputs: Torque [N/m], pitch angle beta [rad], wind speed V [m/s]
% BwithBeta = Beval; %zeros(1,size(B,2))];
% BwithBeta(:,1) = 0;
% BwithBeta = [BwithBeta; zeros(1,size(Beval,2))]; % For new state V
% BwithBeta(end,1) = 1/tau; % For new state V

% add v as state
% AwithV = [AwithBeta, BwithBeta(:,1); [zeros(1,size(AwithBeta,2)),-1/tau_g]];
% BwithV = BwithBeta; %zeros(1,size(B,2))];
% BwithV(:,1) = 0;
% BwithV = [BwithV; zeros(1,size(BwithV,2))]; % For new state V
% BwithV(end,1) = 1/tau_g; % For new state V
AwithV = Aeval;
BwithV = Beval;

%% Calculate discrete forms of A,B
A2 = AwithV*AwithV;
A3 = AwithV*A2;
Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*AwithV+DT*eye(size(AwithV));


%% Finish calculation of A,B
A = Apre*AwithV + eye(size(AwithV));
B = Apre*BwithV;


end
