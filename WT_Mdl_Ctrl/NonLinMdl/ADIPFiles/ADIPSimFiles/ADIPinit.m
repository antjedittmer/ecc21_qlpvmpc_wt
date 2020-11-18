%% Initialize physical and system constants

% Pablo S.G. Cisneros, Herber Werner, ICS TUHH

% initial angles and angle rates
x0 = [0;0;0;0];
gainDisturbance = 0;

% Switch disturbance acting on outputs included
if ~(exist('dOut','var') == 1)
    dOut = 0;
end

if dOut
    x_ke = [x0; 0 0]; %#ok<UNRCH>
    P_ke = zeros(6);
else
    x_ke = [x0; 0];
    P_ke = zeros(5);
end


% Tuning parameter qLPV observer and controller
Qe = 1000; % process noise covariance: 1000 * eye(5); 
Re = 0.001; % measurement noise covariance

N = 40;
Q = [200,1000,0.1,10,0]; %[1000 1000 0.1 10 0];
R = 1000; %kron(eye(N),2000);
Xko_init = kron(ones(N,1),x_ke);

% gravitational constant
g = 9.81;

% ADIP physical parameters
m1 = 0.115; % Arm mass
mh = 0.13; % Encoder 2 mass
m2 = 0.0731; % Pendulum mass
M = mh+m2; 

l1 = .158; % Arm length 
lg1 = l1/2; % Arm length (center of gravity)

l2 = .078;  % pendulum length
I2 = 1.494*m2*l2^2;

fv1 = 1.3e-3; % Arm friction coefficient
fv2 = 2.2e-5; % Arm friction coefficient

%% Initialize sysConstants for usage in constant block
%sysConst = [l1,l2,m1,m2,g,fv1,fv2];
sysConst = [...
l1 ; % = 5.875 *.0254;
l2 ; % =6.8*.0254/1.5; %equivalent mass position
m1 ; % ==.19;
m2 ; % == 0.08;
g ; % == 9.81;
fv1 ; % ==4.4e-3; % fv1=3.4e-3;
fv2 ]; %= 3.4e-5;