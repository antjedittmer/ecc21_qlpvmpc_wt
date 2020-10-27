%clear; clc; close all

%% Load lambda, beta, torque coefficient and assign parameters
load 'NREL5MW_linearised'
load 'NREL5MW_CPdata_new3'

% Size states, inputs outputs
nx = 4; % number states plant: ThetaS,omegaR,omegaG,Tg
ny = 2 ; % number of outputs plant omegaG,Tg
nxA = nx +ny; %number states augmented plant 
nu = 2; % plant inputs: Tg_ref, beta_ref

% Initialize WECS parameter 
rb = 63;
Ng = 97;

% Initialize simulation time
DT = 0.008;
Aj = nan(nxA, nxA, length(Lin.V));
Bj = nan(nxA, nu, length(Lin.V));
C = [zeros(ny,nx-ny), diag([Ng, ones(1,ny-1)])];

% Create matrices for LPV ss array
for index = 1:length(Lin.V) 
    
    %% Variables dependent on wind speed V
    Vbar = Lin.V(index);
    Lambda_bar = Lin.RSpeed(index)*rb/Vbar;
    omegabar = Lin.RSpeed(index);
    beta_bar = Lin.Pitch(index);

    %% Grid points
    rho = [omegabar,beta_bar,Vbar];
    [Aj(:,:,index),Bj(:,:,index)] = ABqLd_WT(DT,rho,Rotor_Lamda,Rotor_Pitch,Rotor_cQ,C);  

end

%% Create ss aray
nxA = nx +ny; %number states augmented plant 
Carray = eye(nxA);
Carray1 = [eye(ny),zeros(ny,nx)];
Darray = zeros(nxA,nu);
Darray1 = zeros(ny,nu);
ssLPV = ss(Aj,Bj,Carray1,Darray1,DT);
impulse(ssLPV(:,:,1))

