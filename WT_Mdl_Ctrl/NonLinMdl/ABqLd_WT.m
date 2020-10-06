function [A,B, Arho, Brho] = ABqLd_WT(DT,rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,C)
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
nx = size(C,2); % Number of states: 4 or 5 (with beta)
[Arho,Brho] = ABqL_WT(rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,nx);
%[Arho,Brho] = ABqL_WT_KMC(rho,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,nx);

%% Calculate new matrices
ny = size(C,1);
nu = size(Brho,2);

% [ydot, xdotdot] = A * [y, xdot] + B * udot
A = [zeros(ny,ny) C; zeros(nx, ny) Arho];
B = [zeros(ny,nu); Brho];

% Discretize matrices
AwithV = A;
BwithV = B;

%% Calculate discrete forms of A,B
A2 = AwithV*AwithV;
A3 = AwithV*A2;
Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*AwithV+DT*eye(size(AwithV));

%% Finish calculation of A,B
A = Apre*AwithV + eye(size(AwithV));
B = Apre*BwithV*1/DT;
 
% Carray1 = [eye(ny),zeros(ny,nx)];
% Darray1 = zeros(ny,nu);
% SYSD = c2d(ss(AwithV,BwithV,Carray1,Darray1),DT,'zoh');
% figure;
% impulse(SYSD); hold on;
% 
% SYSD1 = ss(A,B,Carray1,Darray1,DT);
% impulse(SYSD1,'r--');

end