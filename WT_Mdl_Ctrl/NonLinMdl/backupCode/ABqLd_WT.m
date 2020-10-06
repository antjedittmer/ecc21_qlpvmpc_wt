function [A,B, Arho, Brho] = ABqLd_WT(Ts,rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,C)
% ABqLd calculates system (A) and input (B) matrix for quasiLPV
% representation from ABgL and augments the states with integrating
% disturbances.
%
% Inputs: * Ts: Sampling time
%         * rho: Scheduling parameter vector: rho = V
%
% Outputs: * A: Augmented system matrix A(rho) at current rho
%          * B: Augmented Input matrix B(rho) at current rho

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%% Call ABqL
[Arho,Brho] = ABqL_WT(Ts,rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ);

%% Calculate new matrices
nx = size(Arho,1);
ny = size(C,1);
nu = size(Brho,2);

% [ydot, xdotdot] = A * [y, xdot] + B * udot
A = [zeros(ny,ny) C; zeros(nx, ny) Arho];
B = [zeros(ny,nu); Brho];
end