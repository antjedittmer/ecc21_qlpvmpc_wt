% Author      : Bindu Sharan
% Date        : 13.07.2021
% Description : Calculate optimal controller without lqr command.
%               Check with lqr 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all; clear

%% Define system and input matrix
mu = -5;
lambda = -10;
A = [mu      0       0; 0 lambda -lambda; 0      0    2*mu];
B = [1 1 1]';

%% Set weighting matrices, calculate weighting matrices and lqr solution
r = 1;
Q =[1 0 0;0 1 0;0 0 0]; % weightage given to states of nonlinear system only.
F1 = lqr(A,B,Q,r);

%% Hamiltonion matrix -lambda to make it stabizable
H = [A, -B* r^(-1) * B';
    -Q, -A']; 

%% Divide into stable and unstable eigenvalues
[V,D] = eig(H);

idxS = diag(D) < 0; % get stable eigenvalues
Vs = V(:,idxS);% get eigenvectors of stable eigenvalues
Vus = V(:,~idxS); % not stable

n11 = sum(idxS);
P = Vs(n11+1:end,:)/(Vs(1:n11,:));
F = r^(-1) * B' * P;

disp('Solution lqr');
disp(F1);
disp('Solution hamiltonian');
disp(F);

