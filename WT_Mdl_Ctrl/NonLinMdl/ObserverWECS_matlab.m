function [xk, Pk] = ObserverWECS_matlab(xdotk_1,theta1234,udot_1,Pk_1,Rotor_Lambda,Rotor_Pitch,Rotor_cQ)
% Observer calculates a quasiLPV Kalman Filter (qLPV KF), similar to an 
% Extended Kalman Filter (EKF)
% 
% The quasiLPV KF calculates:
% * quasiLPV state estimate (EKF: non-linear model)
% * quasiLPV model covariances (EKF: linearized model)
%
% Inputs: * xk_1: quasiLPV KF state estimates, delayed by one sample
%         * theta1234: omega_g [rad/s] : 3rd state, first output
%         * theta2: T_g [Nm]: 4th state second output
%         * udot_1: quasiLPVMPC controller output gradient, delayed by one sample
%         * Pk_1: Scheduling trajectory, delayed by one sample. the
%         scheduling trajectory contains [y, xdot] with y = [omega_g,T_g]
%         and xdot = [Thetadot_s, Omegadot_r,Omegadot_g,Tdot_g]
%         
% Outputs: * xk: quasiLPV KF state estimates
%          * Pk: Scheduling trajectory [y, xdot]
%          * Vest based on Vest_old + udot = 

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%#codegen
%% Define noise covariance matrices and set filter for GenSpd and GenTq
% States: Theta_s [rad], omega_r [rad/s], omega_g [rad/s], T_g [Nm]
% Input: V [m/s], beta [deg], T_g_ref [Nm]
% Output: omega_g [rad/s], T_g [Nm]
DT = 0.0080;
% theta1 = theta12(1);
% theta2 = theta12(2);
nxA = 10; %number states augmented plant
nx = 6; % number states plant
ny = 4; % number of outputs plant
PkEye =  eye(nxA);
Ng = 97;

Qe = 1000 * PkEye; % process noise covariance: 1000 * eye(5); 
Re = 0.1; % measurement noise covariance
% C = [1 0 0 0 0; 0 1 0 0 0]; % get angles theta1 and theta2, dismiss rates
C1 = [zeros(ny,nx-ny), diag([Ng; ones(ny-1,1)])]; % get GenSpeed, GenTq, V
rho = xdotk_1([1,3,4]); % Tg not included: [Omega_G [1/s], beta[rad], V [m/s]];

%% States for model
% xA = [y, xdot] with
% y = [omegaG, Tg, beta, V]
% x = [ThetaS, OmegaR, OmegaG, Tg, beta, V]
% xdot = [ThetaS, OmegaR, OmegaG, Tg, beta, V]';
% u = [V,beta_ref,TgRef]; Input: V [m/s], beta [rad], T_g_ref [Nm] (in HSS)

%% Calculate quasiLPV model with current state x = [theta', thetad']
% %[A,B] = ABqLd(Ts,rho,sysConst,act) with rho = [theta1,theta1d]
cQ = Rotor_cQ(:,1:size(Rotor_Lambda,1)); 
[A,B] = ABqLd_WT(DT,rho,Rotor_Lambda,Rotor_Pitch,cQ,C1);
C = [eye(4),zeros(4,6)];

%% Prediction step quasiLPV KF
xdodotk_pre = A*xdotk_1+B*udot_1; %ToDo
Pk_pre =  A*Pk_1*A'+ Qe * eye(size(A,1));

%% Update step quasiLPV KF (see Paper section 3.1 Observer Design)
zk = theta1234'; % zk plant output: GenSpeed, TauG, beta, V
yk = zk - C*xdodotk_pre; % Paper: ytilda_k = y_k - C xhead_{k|k-1} 
Sk = C*Pk_pre*C' + Re * eye(size(C,1)); % Paper: Sk = (C P_{k|k-1} C^T +R_e)
Kk = Pk_pre*C'/Sk; % Paper: K_k = P_{k|k-1}C^T Sk^-1

Pk = (PkEye  - Kk*C)*Pk_pre; % Paper: P_{k|k} = (I - K_k C) P_{k|k-1}

%% Limit scheduling parameters and changes in scheduling parameters
xk = xdodotk_pre + Kk*yk; % Paper: xhead_{k|k} = xhead_{k|k+1} + K_k ytilda_k

% Calculate udot based on xk
%uDot =-pinv(Brho)/ Arho * xk(3:6); %-Brho\ Arho * xk(3:6);
% vDot = - (Arho(2,:)*xk(3:6))/Brho(2,1)
% Vest = Vest + vDot;

