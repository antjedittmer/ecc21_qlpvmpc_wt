function [u,Xk,U,mu,us,tictoc]= qLPVMPC_WECS_matlab(r,xk,Xk,mu_old,lambdaVec,pitchVec,cQ,schedParam,U,DT)
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
N = 20;
q_ = [10 1 1 1 10 1 0]; %q_ = [100 1000 1 0.01 0]; 7.351e+02 1.000e+00 1.938e+09 1.625e+08 1.491e+08 1.657e+03
r_ = [1000,1];
R_ = kron(eye(N),diag(r_)); %R_ = 1000*eye(N); 
p = 100*q_; %P = 1000*diag(q_);

Ulim = repmat([8 * pi/180, 20000]',2*N,1); %rate constraints 

Control_Pitch_Max = 90 /180 * pi;
Control_Pitch_Min = -2 /180 * pi;
Control_Torque_Limit = 47403;
Control_Torque_Min = 200;
schedParam(2) = min(max(schedParam(2), Control_Pitch_Min),Control_Pitch_Max);

% isBetaState = 0;
% nx = 4 + isBetaState; % number states plant: ThetaS,omegaR,omegaG,Tg
% ny = 2 + isBetaState; % number of outputs plant omegaG,Tg
nxA = 7;% nx + ny + 1; %number states augmented plant 
nu = 2; % plant inputs: Tg_ref, beta_ref

%% Shift Xk
Xk(1:end-nxA) = Xk(nxA+1:end); 
%Xk(1:nxA)

%% Calculate input u_dot_s = 0  and steady state x = r
us = [0;0]; %xuss(nxA+1:nxA+2); % xuss = [eye(size(Ass))-Ass -Bss;[1 0 0 0] 0] \ [Bss*d;reference];

iterations = 1; % omega_ref(Vest), Tg_ref(Vest)
refVec = zeros(nxA,1);
refVec(1:2) = r;
ref = repmat(refVec,N,1); %ref = kron(ones(N,1),[r; 0; 0; 0;0]);

for i = 1:iterations
    [L,S,H,g] = HSqLd_WT(N,q_,R_,p,xk,Xk,DT,ref,us,lambdaVec,pitchVec,cQ,schedParam,U);
    
    %% %%%%%%Turn QP into LCP and solve it%%%%%%%%%
    Ac_ = [-eye(N*nu);eye(N*nu)];
    bc_ = - Ulim;
    Hi = eye(size(H,1))/H;
    M = Ac_*Hi*Ac_';
    q = -Ac_*Hi*g - bc_;
 
    mu = lemke(M,q,mu_old);
 
    U = Hi*(Ac_'*mu-g);
    U(isnan(U)) = 0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Xk = L*xk+S*U;
end

u =  [min( max(U(1) + schedParam(2), Control_Pitch_Min),Control_Pitch_Max),... % pitch
min(max(U(2) + Xk(2), Control_Torque_Min),Control_Torque_Limit)]'; % torque

%%%%%%%%%%% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h2(Uk) < 0, at current estimate xl
% <=> min( x'H(xl)x + g(xl)x s.t. dh1(xl)x+h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with 
% Jk = (Lx+SUk-Xss)'Q(Lx+SUk-Xss) +(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss) 

diffX = Xk(1:(end - nxA)) - ref(1:(end -nxA));
diffXMat = sum(reshape(diffX.^2, nxA, N-1),2);

UMat = sum(reshape(U.^2, 2, N),2);

diffXEnd = Xk((end - nxA +1):end) - ref((end - nxA +1):end);

% Q_ = kron(eye(N-1),diag(q_));
% Qcomp = diffX' * Q_ * diffX;
% Pcomp = diffXEnd' * diag(p)* diffXEnd;
% Rcomp = U' * R_ * U;

Qcomp = diffXMat' .* q_;
Pcomp = diffXEnd' .* p;
Rcomp = UMat' .* r_;


tictoc = toc;