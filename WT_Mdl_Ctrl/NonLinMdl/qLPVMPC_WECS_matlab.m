function [u,Xko,Udot,mu,us,tictoc]= qLPVMPC_WECS_matlab(r,xk,Xk,mu_old,lambdaVec,pitchVec,cQ,schedParam,Udot,DT)
% MATLAB_Function calculates controller output u. It turns the QP into a 
% linear complementarity problem (LCP) and solves it with Lemke's algorithm.
%
% Inputs: * r: reference
%         * xk: quasiLPV KF state estimates
%         * Xk: Xko(= L*xk+S*U) delayed by one sample timestep
%         * mu_old: mu delayed by one sample timestep
%         * lambdaVec: lambda vector for cQ LUT 
%         * pitchVec: pitch vector for cQ LUT 
%         * cQ: torque coefficients LUT matrix
%         * schedParam: current values [omegaR, beta, V]
%         * U: U (from Xko = L*xk+S*U) delayed by one sample timestep
%
% Outputs: * u: controller output at current time step
%          * Xko: Future state trajectory (L*xk+S*U)
%          * U: Future input trajectory 
%          * mu: Result Lemke's algorithm/LCP: M mu +q >=0, mu>=0, mu'(Mmu+q)=0
%          * us: xuss(end) from xuss = [eye(4)-Ass -Bss;[1 0 0 0] 0]\[Bss*d;1]
%          * tictoc: tic-toc result

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%#codegen
coder.extrinsic('tic');
coder.extrinsic('toc');
tictoc = 0;
tic

%% Init N, Q, R
N = 40;
%q_ = [100 100 zeros(1,4)]; %q_ = [100 1000 1 0.01 0];
q_ = [100 100 zeros(1,2) 1000 100]; %q_ = [100 1000 1 0.01 0];
R_ = kron(eye(N),diag([1000,10])); %R_ = 1000*eye(N); 
p = 1000*q_; %P = 1000*diag(q_);

Ulim = repmat([8 * pi/180, 20000]',2*N,1); %rate constraints 

Control_Pitch_Max = 90 /180 * pi;
Control_Pitch_Min = -2 /180 * pi;
Control_Torque_Limit = 47403;
Control_Torque_Min = 200;

% number states plant: ThetaS,omegaR,omegaG,Tg, number of outputs plant omegaG,Tg
nxA = length(xk); %number states augmented plant 
nu = 2; % plant inputs: Tg_ref, beta_ref
nx = 4; % number states plant: ThetaS,omegaR,omegaG,Tg
ny = nxA - nx; % number of outputs plant omegaG,Tg

%% Shift Xk and shift Udot
Xk(1:end-nxA) = Xk(nxA+1:end); 
Udot(1:end-nu) = Udot(nu+1:end); 

%% Create augmented system: C matrix orig. mdl, A, B matrices augmented
Ng = 97;
C1 = [zeros(ny,nx-ny), diag([Ng, ones(1,ny-1)])];
us = zeros(nu,1); % xuss = [eye(size(Ass))-Ass -Bss;[1 0 0 0] 0] \ [Bss*d;reference];

iterations = 1; % omega_ref(Vest), Tg_ref(Vest)
refVec = zeros(nxA,1);
refVec(1:2) = r;
ref = repmat(refVec,N,1); %ref = kron(ones(N,1),[r; 0; 0; 0;0]);

for i = 1:iterations
    [L,S,H,g] = HSqLd_WT(N,q_,R_,p,xk,Xk,DT,ref,us,lambdaVec,pitchVec,cQ,C1,schedParam,Udot);
    
    %% %%%%%%Turn QP into LCP and solve it%%%%%%%%%
    Ac_ = [-eye(N*nu);eye(N*nu)];
    bc_ = - Ulim;
    Hi = eye(size(H,1))/H;
    M = Ac_*Hi*Ac_';
    q = -Ac_*Hi*g - bc_;
 
    mu = lemke(M,q,mu_old);
 
    Udot = Hi*(Ac_'*mu-g);
    Udot(isnan(Udot)) = 0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Xko = L*xk+S*Udot;
end

u =  [min(max(Udot(1) + schedParam(2), Control_Pitch_Min),Control_Pitch_Max),... % pitch
min(max(Udot(2) + Xk(2), Control_Torque_Min),Control_Torque_Limit)]; % torque

tictoc = toc;