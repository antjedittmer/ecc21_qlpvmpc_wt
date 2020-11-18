function [A,B] = ABqL(Ts,rho)
% ABqL calculates system (A) and input (B) matrix for quasiLPV
% representation
% 
% Inputs: * Ts: Sampling time
%         * rho: Scheduling parameter vector: rho = [theta1,theta1d]
%
% Outputs: * A: System matrix A(rho) at current rho, discretised with Ts
%          * B: Input matrix B(rho) at current rho, discretised with Ts

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

%% Parameters
l1 = 5.875 *.0254;
l2 = 6.8*.0254/1.5; %equivalent mass position

m1 = .19;
m2 = 0.08;
g = 9.81;

fv1 =4.4e-3; % fv1=3.4e-3;
fv2 = 3.4e-5;

%% cos and sin terms for current angles
cost = cos(rho(1));
sinttd = sin(rho(1))*rho(2);
sinct = sinc(rho(1)/pi);

%% M, D, K matrices
M = [(m1+m2)*l1^2        m2*l1*l2*cost;...
    m2*l1*l2*cost  m2*l2^2];

D = [fv1 0;...
    -m2*l1*l2*sinttd fv2];

K = [-(m1+m2)*l1*g*sinct 0;...
    0                       -m2*l2*g];

%% Calculate Aeval from M,D,K
act = 1;
MaskB = [eye(act);zeros(act)];
Mi = inv(M);

Aeval = [zeros(2) eye(2);...
    -Mi*K -Mi*D];

%% Calculate Beval from M,D,K
Beval = [zeros(2,2);Mi]*MaskB;
A2 = Aeval*Aeval;
A3 = Aeval*A2;
%II = eye(4);
Apre = (1.0/24.0)*Ts^4*A3+(1.0/6.0)*Ts^3*A2+0.5*Ts^2*Aeval+Ts*eye(4);

%% Calculate A,B
A = Apre*Aeval + eye(4);
B = Apre*Beval;

end
