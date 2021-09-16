%% Proof-of concept for unknown input observer for WECS
close all;
clear;
addpath(genpath(pwd))
parentDir = fileparts(pwd);

%% Define user input
rhoP = 18; %Wind speed: Scheduling parameter, alternative [V,omega,beta]
useFASTmdl = 0; % 0 (Default) 5th order model used, 1 9th order
useFASTwind = 1; % 0 Wind speed step function, 1: NTW Wind from FASTTool
varWind = 0; % Wind variance
multTau = 0.0175; %Factor for time constant wind speed low pass ( tau = 1/DT)
multVar = 1; % Change variance NTW

%% Set-up wind test signal
% Get time and input vectors to simulate with lsim
DT = 0.008; % sample time
t = 0:DT:100; % time vector

if useFASTwind == 0 %create step function
    t1L = round(length(t)/2);
    V = [rhoP*ones(1,t1L),(rhoP+2) * ones(1,length(t) - t1L)]';
else % cut wind to length of time vector
    dataInPath = fullfile(parentDir ,'dataIn');
    load(fullfile(dataInPath,'OutDataWind18NTW.mat'),'OutTable');
    
    windTest = OutTable.Wind1VelX;
    meanWindTest = mean(windTest); % varWindTest = var(windTest);
    OutTable.Wind1VelX = multVar*(windTest - meanWindTest) + meanWindTest;
    
    V = OutTable.Wind1VelX(1:length(t));
    
end

d1 = varWind*randn(size(V));
d = V + d1;  % wind input load('d'); d = d';

%% Get parameters from Data Dictionary
% The data dictionaries of  WECS test harness with the lookup table (LUT)
LinPointsMFile;
MdlMFile;

%% Set up control inputs as feedforward  and initial state rotor speed
% u = [Tg_ref beta_ref V]'

% Get torque, pitch and and rotor speed corresponding to wind input
u1 = nan(size(d));
u2 = nan(size(d));
rhoP1 = nan(size(d));
for idx = 1: length(d)
    idxVec = abs(Lin_points.V - d(idx)); % difference to wind speed
    idxNearestV =  min(idxVec) == idxVec; % index of nearest wind speed
    u1(idx) = Lin_points.Torque(idxNearestV )*10^3; % Torque [Nm]
    u2(idx) = Lin_points.Pitch(idxNearestV ); % Pitch [rad]
    rhoP1(idx) = Lin_points.RSpeed(idxNearestV ); %Rotor speed [rad/s]
end

%% Set up Simulation model
% The input and output vectors are:
%  u = [Tg_ref beta_ref V]'
% * sysOut.InputName' {'Torque [kNm]'}    {'\beta [rad]'}    {'V [m/s]'}
% y = [yt_dot zeta_dot omega_r omega_gr]'
% {'omega [rad/s]'}    {'ydotdot_{fa} [m/s^2]'}    {'beta'}    {'Tg'}

if useFASTmdl == 1
    % Model form MdlLin folder:
    % Mdl2 is a 9th order model with state vector:
    % x = [yt zeta xt theta_s yt_dot zeta_dot xt_dot omega_r omega_gr]'
    
    % add path to LinMdl folder and run create linear model for windspeed
    addpath(fullfile(parentDir,'LinMdl'));
    sysOut = createLinearModels(rhoP-3,1);
    
    % State names (input and output names already created)
    sysOut.StateName = {'yt [m]', 'zeta [m]', 'xt [m]', 'theta_s [rad]', ...
        'yt_dot [m/s]', 'zeta_dot [m/s]', 'xt_dot [m/s]',...
        'omega_r [rad/s]', 'omega_gr [rad/s]','Tg [kNm]','beta [rad]'};
    
    % Assign matrices
    Av = sysOut.A; Bu = sysOut.B(:,1:2); Cv = sysOut.C(1:4,:); Ev = sysOut.B(:,3);
    
    % Initial steady state and output (testX = -inv(Av) * Bv * u0'; testX(1:3)
    x0 = [zeros(1,7),rhoP1(1),rhoP1(1),u1(1),u2(1),d(1)]';
    y0 = x0(8:end-1); %{'omega_r','omega_gr','Tg','beta','V'};
    
else 
    
    % System and input matrix
    [Av,Bv] = ABqL_WT_beta(d(1),lambda,beta,Cq);
    
    %idxState = 1:3;   
    Av = Av(:,:); 
    Bu = Bv(:,2:3);
    nx = size(Av,2);
    Cv = [zeros(nx-1,1),eye(nx-1)]; %
    Ev = Bv(:,1);
    
    % Initial steady state and output (testX = -inv(Av) * Bv * u0'; testX(1:3)
    x0 = [0,rhoP1(1),rhoP1(1),u1(1),u2(1),d(1)]';
    
    y0 = x0(2:end-1);
    
    % State and output names
    sysOut.StateName = {'theta_s [rad]','omega_r [rad/s]',...
        'omega_{gr} [rad/s]','Tg [Nm]','beta [rad]'};
    sysOut.OutputName = [sysOut.StateName(2:5),{'V [m/s]'}];
    
%     sysTest = ss(Av,Bv,C,zeros(size(C,1),size(Bv,2)));
%     u0 = [V(1),u1(1),u2(1)];
%     Atest = - Bv*u0' / x0(1:5);
    
end

%% Check rank
fprintf('WECS: rank(C*E) is %d; rank(E) is %d\n', rank(Cv*Ev),rank(Ev));
nx = size(Av,1);
ny = size(Cv,1);

% Also create an additional state for the disturbance
tau1 = 1/DT* multTau;
A =[[Av, zeros(nx,1)]; [zeros(1,nx),-tau1]];
C = [Cv, zeros(ny,1)]; %[zeros(1,nx),1]];
E = [Ev;tau1]; %xd_dot = -xd + Ed
B = [Bu; zeros(1,2)];
fprintf('rank(C*E) is %d; rank(E) is %d\n', rank(C*E),rank(E));

%% Calculate observer matrices
H = E*inv((C*E)'*(C*E))*(C*E)' %#ok<MINV,*NOPTS> display output, inv is fast enough
nx = size(A,1);
ny = size(C,1);
T = eye(nx)-H*C
A1 = T*A

% Pole placement with Q and R
Qe = B*B'; %
Re = C*C';
%K1 = lqr(A',C',Qe,Re)'; % lqr not supported in Simulink

Ham = [A', -C'* Re^(-1) * C; % Hamiltonion matrix -lambda to make it stabizable
    -Qe, -A]; 
[V,D] = eig(Ham); % Divide into stable and unstable eigenvalues
idxS = diag(D) < 0; % get stable eigenvalues
Vs = V(:,idxS);% get eigenvectors of stable eigenvalues
Vus = V(:,~idxS); % not stable
n11 = sum(idxS);
P = Vs(n11+1:end,:)/(Vs(1:n11,:));
K1 = real(Re^(-1) * C * P)';

% Calculate F and K
F = A1 - K1*C; %T*A -K1*C
K = K1 + F*H;
T = eye(size(H*C)) - H*C

% Sys with no ctrl input: xdot = Ax + Ed, y = Cx
nd = size(E,2); %nu = 0;
sys = c2d(ss(A,[E,B],C,zeros(ny,nd+2)),DT,'zoh'); % ss(A,[E,B],C,zeros(ny,nd+2));

% Est with no control input: zdot = Fz + Ky, xhat = z + Hy
sysEst = c2d(ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]),DT,'zoh'); %ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]); %
sysEstCont = ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]);

in1 = [d,u1,u2];
diffIn = diff(in1);

in0 = in1 - in1(1,:); %For test
[yhat,~,xhat] = lsim(sys,in0,t);
y1 = yhat + y0';
x1 = xhat + x0';

[dy,~,dx] = lsim(sys,diffIn,t(1:end-1));
ytemp = y0'+cumsum(dy);
y = [ytemp;ytemp(end,:)];
xtemp = x0' + cumsum(dx);
x = [xtemp;xtemp(end,:)];

inEst = [u1,u2,y];
inEst0 = inEst - inEst(1,:);
z0 = T*x0*0; %(eye(size(H*C))- H*C);
[xesthat,~,zhat] = lsim(sysEst,inEst0,t,z0);
z1 = zhat + z0';
xest1 = xesthat + x0';

inEst = [u1,u2,y];
diffInEnst = diff(inEst);
[dXest,~,z] = lsim(sysEst,diffInEnst,t(1:end-1));
xtempEst = x0' + cumsum(dXest);
xest = [xtempEst;xtempEst(end,:)];

figure(1)
cl = colormap('lines');

ny = size(y,2);
for idx = 1:ny
    subplot(ny, 1,idx)
    plot(t,y(:,idx), 'Color',cl(1,:),'Linewidth',1);
    hold on;
     plot(t,y1(:,idx),'Color',cl(2,:), 'Linewidth',1,'Linestyle','--');
    if idx == 1
        title('WECS outputs');
    end
    axis tight;
    grid on;
    ylabel(sysOut.OutputName{idx});
end
xlabel('time [s]');


figure(2)
idxEnd = size(x,2)-1;
for idx = 1:idxEnd
    subplot(idxEnd,1,idx)
    plot(t,x(:,idx),'-','Color',cl(1,:));
    hold on;
    plot(t,xest(:,idx),'k--');
    axis tight; grid on;
    if idx == 1
        title('WECS states and estimates');
    end
    axis tight;
    grid on;
    ylabel(sysOut.StateName{idx});
end
xlabel('time [s]');

figure(3)
idx = size(x,2);

plot(t, d,'k','Linewidth',1);
hold on;
plot(t,x(:,idx),'-','Color',cl(1,:),'Linewidth',1);
plot(t,xest(:,idx),'--','Color',cl(2,:),'Linewidth',1);
legend('wind','wind state','wind est');
ylabel('wind [m/s]')
xlabel('time [s]');
title('Comparison Measured and Estimated Wind');
axis tight;
grid on;
print(gcf,fullfile('figDir','windEst'), '-dpng');

vafState = 1 - var(x(:,end)-d)/var(d)
vafStateEst = 1 - var(xest(:,end)-d)/var(d)
