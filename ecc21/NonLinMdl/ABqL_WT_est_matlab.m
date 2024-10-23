function [Ak,Bk,Ck,Dk]  = ABqL_WT_est_matlab(rhoP,lambda,beta,Cq,DT)
% ABqL calculates system (A) and input (B) matrix for quasiLPV
% representation for ve
%
% Inputs: * Ts: Sampling time
%         * AwithV: Scheduling parameter vector: AwithV = V
%
% Outputs: * A: System matrix A(AwithV) at current estimated air speed AwithV, discretised with Ts
%          * B: Input matrix B(AwithV) at current AwithV, discretised with Ts

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for WECS simulation: Antje Dittmer

persistent lambdaVec pitchVec CqVec Cqdlambda Cqdbeta

if isempty(lambdaVec)
    lambdaVec =  lambda;
    pitchVec = beta;
    CqVec = Cq;
    %% Compute gradients %%[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
    [Cqdlambda,Cqdbeta] = gradient(CqVec,mean(diff(lambda)),mean(diff(beta))); %dcQ/dbetaDeg * dbetaDeg/dbetaRad *betaRAd
end

%% Parameters
if length(rhoP) == 1
    Vbar = rhoP; % Vbar
    LUT_V = 4:25; % LUT table rows
    
    LUT_omega = ...%[0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)];
        [0.7458832277625318 0.77759738173519954 0.820732502042672 ...
        0.87379091011167287 0.93538972625666472 1.0429417531956238 ...
        1.1591854391157337 1.2371015764176956 1.2554479267196832 ...
        1.2547857982936621 1.2545286368009017 1.2546177195393646 ...
        1.2545376249695805 1.2545426663335242 1.2545426663335242 ...
        1.2545426663335242 1.2545426663335242 1.2545426663335242 ...
        1.2545426663335242 1.2545426663335242 1.2545426663335242 ...
        1.2545426663335242];
    
    LUT_beta =  [0 0 0 0 0 0 0 0 0.064193209888351435 0.11259119004615419 ...
        0.14929546421559495 0.18081611050661253 0.20926497731412011 ...
        0.23701571242082992 0.26197392072434883 0.28536133270107289 ...
        0.30770154712660031 0.32934362985133 0.35028758087526196 ...
        0.37035886727319672 0.38973202197033374 0.40858157789187255 ...
        ];
    
    
    % Get omegabar and beta
    omegabar = interp1(LUT_V,LUT_omega,Vbar,'pchip');
    beta_bar = interp1(LUT_V,LUT_beta,Vbar,'pchip');
    
elseif length(rhoP) == 3
    % rho = xdotk_1([1,3,4]); % Tg not included: [Omega [1/s], beta[rad], V [m/s]];
    Vbar = rhoP(3);
    omegabar = rhoP(1); % pass OmegaR OmegaG/Ng as approximation of omegaR
    beta_bar = rhoP(2);
    
else
    error('Unexpected length of scheduling parameter rho')
end

%% Physical and model constants
rho = 1.2250;
rb = 63; % m blade radius
Ng = 97; %Gearbox ratio
Jr = (115926 + 3 * 11.776e6) *1; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
Jgr = Jg * Ng^2; % Inertia of the generator in LSS
Ks = 867637000; %Nm/rad Stiffness of the transmission
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
tau_g = 0.01; % for pitch actuator and airspeed

%% Lambda
Lambda_bar = omegabar*rb/Vbar;
dlambdadomega= rb/Vbar;
dlambdadV = - omegabar*rb/Vbar^2;

%% Grid points
Cqb = interp2(lambdaVec,pitchVec,Cq,Lambda_bar,beta_bar);
dCqdlambdab = interp2(lambdaVec,pitchVec,Cqdlambda,Lambda_bar,beta_bar);
dCqdbetab = interp2(lambdaVec,pitchVec,Cqdbeta,Lambda_bar,beta_bar);
%% Tr=0.5*pi*rb^3*Cq*v^2
dTrdomega = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;% br = dTr/domega
dTrdbeta = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;% Krb = dTr/dbeta
dTrdV =  2*0.5*rho*pi*rb^3*Vbar*Cqb + ...
    0.5*rho*pi*rb^3*Vbar^2*dCqdlambdab*dlambdadV;

%% Create system matrix A:
% States: [Theta_s [rad], omega_r [rad/s], omega_gr [rad/s], T_g [Nm]]
% States for model xA = [y, xdot] with y = [Tg, beta]
% x = [ThetaS, OmegaR, OmegaG, Tg, beta,]
Atemp = [[0,1,-1,0] ;
    [-Ks/Jr, (-Bs + dTrdomega)/Jr, Bs/Jr,0];
    Ks/Jgr, Bs/Jgr, -Bs/Jgr, -Ng/Jgr;
    zeros(1,3), -1/tau_g];

%% Create input matrix B
% Input orig model: V [m/s], T_g_ref [Nm], beta [rad],
% -> vel model: Vdot [m/s^2]  , T_g_refDot [Nm/s]
% for beta: Vdot [m/s^2]  T_g_refDot [Nm/s] betaDot [rad/s],
Btemp = [zeros(1,3); dTrdV/Jr, 0 dTrdbeta/Jr; zeros(1,3);  0,  1/tau_g,0];

% With Beta
AvAll = [Atemp,Btemp(:,3); zeros(1, size(Atemp,2)), -1/tau_g];
BvAll = [Btemp; zeros(1, size(Btemp,2))];
BvAll(:,3) = 0;
BvAll(end,3) = 1/tau_g;


idxS = 3;
Av = AvAll(1:idxS,1:idxS);
Bv = [BvAll(1:idxS,1),AvAll(1:idxS,idxS+1), AvAll(1:idxS,idxS+2)];


%% Check rank
Bu = Bv(:,2:3);
Bd = Bv(:,1);
nxEval = size(Av,2);

% Also create an additional state for the disturbance
tau1 =  1/DT;
A =[[Av, zeros(nxEval,1)]; [zeros(1,nxEval),-tau1]];
E = [Bd;tau1]; %xd_dot = -xd + Ed
B = [Bu; zeros(1,2)];

%% Calculate observer matrices
Ceval = [zeros(nxEval-1,1),eye(nxEval-1)]; %
nyEval = size(Ceval,1);

C = [Ceval, zeros(nyEval,1)]; %[zeros(1,nx),1]];
H = E*inv((C*E)'*(C*E))*(C*E)'; %#ok<MINV,*NOPTS> display output, inv is fast enough
nx = size(A,1);
T = eye(nx)-H*C;
A1 = T*A;

% Pole placement with Q and R
Qe = [E,B]*[E,B]'; %eye(size(A)); %B*B'; %
Re = C*C';
%K1 = lqr(A',C',Qe,Re)'; % lqr not supported in Simulink

Ham = [A', -C'* Re^(-1) * C; % Hamiltonion matrix -lambda to make it stabizable
    -Qe, -A]; 
[V,D] = eig(Ham); % Divide into stable and unstable eigenvalues
idxS = diag(D) < 0; % get stable eigenvalues
Vs = V(:,idxS);% get eigenvectors of stable eigenvalues
n11 = sum(idxS);
P = Vs(n11+1:end,:)/(Vs(1:n11,:));
K1 = real(Re^(-1) * C * P)';

% Calculate F and K
F = A1 - K1*C; %T*A -K1*C
K = K1 + F*H;
T = eye(size(H*C)) - H*C;

% Sys with no ctrl input: xdot = Ax + Ed, y = Cx  Est with no control input: zdot = Fz + Ky, xhat = z + Hy
% sysEst = ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]); %

%% Calculate discrete forms of A,B
A2 = F*F;
A3 = F*A2;
Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*F+DT*eye(size(F));

%% Finish calculation of A,B
Ak = Apre*F + eye(size(F));
Bk = Apre*[T*B,K];
Ck = eye(size(F));
Dk = [zeros(size(T*B)),H];

