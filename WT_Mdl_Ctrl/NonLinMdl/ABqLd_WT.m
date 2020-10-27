function [A,B] = ABqLd_WT(DT,rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,C)
% ABqLd calculates system (A) and input (B) matrix for quasiLPV
% representation from ABgL and augments the states with integrating
% disturbances.
%
% Inputs: * DT: Sampling time
%         * rho: Scheduling parameter vector: rho = V
%
% Outputs: * A: Augmented system matrix A(rho) at current rho
%          * B: Augmented Input matrix B(rho) at current rho

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%% Call ABqL
[Arho,Brho] = ABqL_WT(rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ);

%% Calculate discrete forms of A,B
A2 = Arho*Arho;
A3 = Arho*A2;
Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*Arho+DT*eye(size(Arho));

%% Finish calculation of A,B
Ak = Apre*Arho + eye(size(Arho));
Bk = Apre*Brho;

%% Calculate new matrices augmented system
ny = 2;
% nu = size(Brho,2);

A = [eye(2) diag([97,1])* Ak(3:4,:);
    zeros(size(Ak,1),ny), Ak];
B = [diag([97,1])* Bk(3:4,:); Bk];

% [ydot, xdotdot] = A * [y, xdot] + B * udot
% A = [zeros(ny) C *1/DT; zeros(nx, ny) Arho];
% B = [zeros(ny,nu); Brho];
 
% nx = 4; Carray = [zeros(ny,ny),eye(ny)];
% Darray1 = zeros(ny);
% SYSD = c2d(ss(Arho,Brho,Carray,Darray1),DT,'zoh');
% figure;
% impulse(SYSD); hold on;
% 
% SYSD1 = ss(A,B,Carray1,Darray1,DT);
% impulse(SYSD1,'r--');

end