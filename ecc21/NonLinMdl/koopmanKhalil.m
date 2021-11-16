% Author      : Bindu Sharan
% Date        : 13.07.2021
% Description : This is script is meant to be used to analyze the
%               Stability of the nonlinear system. First the nonlinear system is lifted into finte order higher dimension.
%               then the remaining terms have been considered as perturbations and theory from Khalil (ch 9)has been
%               applied to find the trade of between the control action and area of action in LQR control
%               Example of a nonlinear model that can be represented
%               as a linear system with perturbations one using the Koopman operator framework. The
%               nonlinear model is
%
%
%               dx1 =       mu*x1      + b11*u + d1,
%               dx2 = lambda*(x2-x1^2) + b21*u + d2,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define system
% using the Koopman framework we can represent the system in the following
% form       dx = Ax + Bu with
%
%              | mu        0        0 |        |0|       |1|
%           A =|  0   lambda  -lambda | ,  B = |1|   or  |0|
%              |  0        0     2*mu |        |0|       |0|.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Are the configurations stabilizable
% There are 4 possible combinations of parameters and 2 possible
% configurations this yields 8 possible models.

% For Configuration 1 i.e. B1 (related eigenvalue is lambda)
% 1- (+ +) i.e. mu>0, lambda>0 ---> Not stabilizable.
% 2- (- -) i.e. mu<0, lambda<0 ---> Stabilizable
% 3- (+ -) i.e. mu>0, lambda<0 ---> Not stabilizable.
% 4- (- +) i.e. mu<0, lambda>0 ---> Stabilizable.


%For configuration 2 i.e. B2 (related eigenvalue is mu)
% 5- (+ +) i.e. mu>0, lambda>0 ---> Not stabilizable.
% 6- (- -) i.e. mu<0, lambda<0 ---> Stabilizable
% 7- (+ -) i.e. mu>0, lambda<0 ---> Stabilizable.
% 8- (- +) i.e. mu<0, lambda>0 ---> Not stabilizable.
%%
clc
close all
clear
n      = 3;
%syms mu lambda;

mu = 0.05;
lambda = 1;

B1 = [1 0 0]';
B2 = [0 1 0]';
N = 2; % N represents the case
switch N
    case 1
        mu = 0.05;
        lambda = 1;
        B = B1;
    case 2
        mu = -0.05;
        lambda = -1;
        B=B1;
    case 3
        mu=0.05;
        lambda=-1;
        B=B1;
    case 4
        mu=-0.05;
        lambda=1;
        B=B1;
    case 5
        mu=0.05;
        lambda=1;
        B=B1;
    case 6
        mu=-0.05;
        lambda=-1;
        B=B2;
    case 7
        mu=0.05;
        lambda=-1;
        B=B2;
    case 8
        mu=-0.05;
        lambda=1;
        B=B2;
    otherwise
        disp('you have entered an invalid value')
end

A = [mu      0       0;...
    0 lambda -lambda;...
    0      0    2*mu];


%% Controllability analysis (Both configurations are not controllable).
isControllableC1 = (rank(ctrb(A,B))==n);
%% Assumptions all observables can be measured
eps11 = 2*mu; eps22 = 2*lambda; eps33 = 4*mu;
A1 = A - [2*mu 0 0;0 2*lambda 0; 0 0 4*mu];
epsilon3=4*mu;
A3 = A -epsilon3*[0 0 0;0 0 0;0 0 1];
% syms r
r = 1;
C=eye(3);
Q=[1 0 0;0 1 0;0 0 0]; % weightage given to states of nonlinear system only.
if N==1
     H=[A1 -B* inv(r)*B';-Q -A1']; %
elseif   N==2
    H= [A - B* inv(r)*B';-Q -A'];% hamiltonion matrix -lambda to make it stabizable
elseif N==3
    H= [A3 -B*inv(r)*B';-Q -A3'];
elseif N==4
    H=[A4 -B*inv(r)*B';-Q -A4'];
end

%% eigen vector calculation (if required)
%syms x1 x2 x3 x4 x5 x6
%X=[x1, x2, x3, x4, x5, x6]'
%[V,D]=eig(H);
%lambda1=D(1,1);
%lambda2=D(2,2);
%lambda3=D(5,5);
%(H-lambda1*eye(6))*X;
%V1=[0 1 0 0 1/2 5/11]';
%(H-lambda1*eye(6))*V1;

%%
[V,D] = eig(H);
idxS = diag(D) < 0; % get stable eigenvalues
Vs = V(:,idxS);% get eigenvectors of stable eigenvalues
Vus = V(:,~idxS); % not stable
Vn = [Vs Vus];
%Dn=Vn\H*Vn;  %why not coming
n11 = sum(idxS);
P = Vs(n11,:)* inv(Vs(1:n11,:));
eig(P);
eig(P(2:3,2:3));
sigma = svd(Q+P*B*inv(r)*B'*P);
F = inv(r)*B'*P;
T= Q+P*B*inv(r)*B'*P;
%E=-(P*A+A'*P);
%%
%test
%A=[1 2 3;4 5 6;7 8 9];
%[V,D]=eig(A);
%Vn=[V(:,2) V(:,3) V(:,1)];%simplify
%Dn=A\Vn*Vn; %inv(Vn)*A*Vn
%Dn1=inv(Vn)*A*Vn;

%%

