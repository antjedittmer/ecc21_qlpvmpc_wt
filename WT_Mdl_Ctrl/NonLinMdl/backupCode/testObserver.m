%% get values
%clear; 
close all;
initModel5MWNRELBianchi;
load('workSpace.mat')
load('OutDataSweep.mat');
OutTable = array2table(OutData,'VariableNames', OutList);
idxWind = contains(OutList, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));


%% 
LUT_V = 4:25;
LUT_omega = [0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)];
LUT_beta = [zeros(1,8),0.0642    0.1126    0.1493    0.1808    0.2093    0.2356    0.2602...
    0.2833    0.3054    0.3269    0.3477    0.3676    0.3868    0.4056];
LUT_torque = [5.0870    7.9450   11.4300   15.5300   20.2400   25.5200   31.3500  ,...
    repmat(43.0900,1,15)]; 

Vbar = 4;
rb = 63; % m blade radius
G = 97;
omegabar = interp1(LUT_V,LUT_omega,Vbar);
Lambda_bar = omegabar*rb/Vbar;
beta_bar = interp1(LUT_V,LUT_beta,Vbar)/pi*180;

cQ = Rotor_cQ(:,1:size(Rotor_Lambda,1));
Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
Torque_g_bar = 0.5*rho*pi*rb^3*Cqb*Vbar^2/G; %in HSS

% theta123 = [omega_g_bar, Torque_g_bar, beta_bar/180*pi]; [1/s,Nm,rad]
theta1234 = [omegabar, Torque_g_bar, beta_bar/180*pi];

udot_1 = [0,0]'; % beta_d, Tg_d

xDdot_k1  = zeros(6,1); %xa: [y,xdot]: 
xDdot_k1(1:4)  = theta1234;
% y = [omegaG, T_g, betabar, Vest]
% xDot = [Thetas_dot, omega_r_dot, omega_g_dot, Tg_dot, betaDot, VestDot]
% u = [Vbar,Tg_ref, beta_ref]
P_ke = zeros(6);
%P_ke(:,1:3) = repmat([theta12, Vbar],8,1);

xDdot_Ve = NaN(100+1,10);
xDdot_Ve(1,:) = xDdot_k1;
for i = 1:1000
    [xDdot_k1,P_ke] = ObserverWECS_matlab(xDdot_k1,theta1234,udot_1,P_ke,Rotor_Lambda,Rotor_Pitch,Rotor_cQ);
    xDdot_Ve(i+1,:) = xDdot_k1;
end

xest = xDdot_Ve;
xCell = {'ThetaDot_s [rad]','omegaDot_r [rad]',...
    'omegaDot_g [rad]','TDot_g','BetaDot_r','V_dot'}; %'V_dot'
xCellTable = {'omega_g', 'T_g', 'ThetaDot_s','omegaDot_r',...
    'omegaDot_g','TDot_g'}; %'V_dot'

figure(1); 
subplot(2,1,1); plot(xest(:,1)/97)
hold on; plot(xest(:,3)*180/pi); plot(xest(:,4))
grid on; axis tight
legend('omega_ [rad/s]','beta [deg]','V [m/s]');
subplot(2,1,2); plot(xest(:,2));legend('T_g [Nm]'); 
grid on; axis tight

figure(2);
for idx = 1:3
    subplot(3,1,idx);
    plot(xest(:,idx+4))
    axis tight; grid on;
    ylabel(xCell{idx});
end

figure(3);
for idx = 1:3
    subplot(3,1,idx);
    plot(xest(:,idx+4+3))
    axis tight; grid on;
    ylabel(xCell{idx+3})
end

r(1) = interp1(LUT_V,LUT_omega,Vbar) * Ng; % 
r(2) = interp1(LUT_V,LUT_torque,Vbar) * 1000;

% Out: omegaG, Tg 
% States: ThetaS, omegaR, omegaG, Tg
% In: beta_ref, Tg_ref;
% Scheduling param rho: V, beta_ref, omegaG/Ng (as approximation of omegaR)

%xAug = [y, xdot]
%xAug = [omegaG, Tg,beta, AwithBeta 

xk = [r(1), r(2), zeros(1,4)];
nu = 2;

N = 40;
% Q = [200,1000,0.1,10,0]; %[1000 1000 0.1 10 0];
% R = 1000; %kron(eye(N),2000);
Xko = kron(ones(N,1),xk);
mu_old = zeros(2*N*nu,1);
Vest = Vbar;

% Vbar = rho(3);
% omegabar = rho(1)/Ng; % OmegaG/Ng as approximation of omegaR
% beta_bar = rho(2)/pi*180;
schedParam = [omegabar, beta_bar, Vbar];
[u,Xk,mu,us,tictoc] = qLPVMPC_WECS_matlab(r,xk,Xko,mu_old,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,schedParam);
