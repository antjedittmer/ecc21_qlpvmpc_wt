function [A,B] = ABqLd(Ts,rho)
% ABqLd calculates system (A) and input (B) matrix for quasiLPV
% representation from ABgL and augments the states with integrating
% disturbances.
%
% Inputs: * Ts: Sampling time
%         * rho: Scheduling parameter vector: rho = [theta1,theta1d]
%
% Outputs: * A: Augmented system matrix A(rho) at current rho
%          * B: Augmented Input matrix B(rho) at current rho

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%% Call ABqL
[Arho,Brho] = ABqL(Ts,rho);

%% Calculate new matrices
nx = size(Arho,1);
nd = 1;
B_d = Brho;

A = [Arho B_d; [zeros(nd,nx) eye(nd)]];
B = [Brho; zeros(nd,1)];

end