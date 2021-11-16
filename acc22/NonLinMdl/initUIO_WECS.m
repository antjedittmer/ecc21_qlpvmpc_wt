%function  initUIO_WECS(idxFig)

%% Proof-of concept for unknown input observer for WECS
%if ~nargin

close all; clear; clc;
restoredefaultpath;
%end

pwd1 = fileparts(mfilename('fullpath'));
idxFig = 10;
addpath(genpath(pwd1))
parentDir = fileparts(pwd1);
mainDir = parentDir;

dataInDir = fullfile( mainDir,'dataIn');
addpath(dataInDir);

% Set path to figure directory
figDir = fullfile(parentDir,'figDirUIO');
if ~isfolder(figDir)
    mkdir(figDir)
end


%% Define user input
rhoP = 18; %Wind speed: Scheduling parameter, alternative [V,omega,beta]
useFASTmdl = 0; % 0 (Default) 5th order model used, 1 9th order
useFASTwind = 1; % 0 Wind speed step function, 1: NTW Wind from FASTTool
varWind = 0; % Wind variance
multTau = 1; %0.0175; %Factor for time constant wind speed low pass ( tau = 1/DT)
multVar = 1; % Change variance NTW
idxS = 3;  % 3 states for default model: Thetas, omegar, omegag

%% Set-up wind test signal
% Get time and input vectors to simulate with lsim
DT = 0.008; % sample time
t = 0:DT:100; % default time vector

if useFASTwind == 0 %create step function
    t1L = round(length(t)/2);
    V = [rhoP*ones(1,t1L),(rhoP+2) * ones(1,length(t) - t1L)]';
    
else % cut wind to length of time vector
    dataInPath = fullfile(parentDir ,'dataIn');
    load(fullfile(dataInPath,'OutDataWind18NTWnoise.mat'),'OutTable');
    %load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');
    
    windTest = OutTable.Wind1VelX;
    meanWindTest = mean(windTest); % varWindTest = var(windTest);
    OutTable.Wind1VelX = multVar*(windTest - meanWindTest) + meanWindTest;
    
    %t = OutTable.Time;
    V = OutTable.Wind1VelX(1:length(t));
    
end

d1 = varWind*randn(size(V));
d = V + d1;  % wind input load('d'); d = d';

%% Get parameters from Data Dictionary
% The data dictionaries of  WECS test harness with the lookup table (LUT)
LinPointsMFile;
MdlMFile;
ControllerMFile;
switchNoise = 0;
switchWind = 1;

%% Set up control inputs as feedforward  and initial state rotor speed
% u = [Tg_ref beta_ref V]'

% Get torque, pitch and and rotor speed corresponding to wind input
u1 = interp1(Lin_points.V,Lin_points.Torque,d)*10^3;
u2 = interp1(Lin_points.V,Lin_points.Pitch,d);
rhoP1 = interp1(Lin_points.V,Lin_points.RSpeed,d);

%% Set up Simulation model
% The input and output vectors are:
% u = [Tg_ref beta_ref V]'
% sysOut.InputName' {'Torque [kNm]'}    {'\beta [rad]'}    {'V [m/s]'}
% y = [yt_dot zeta_dot omega_r omega_gr]'
% {'omega [rad/s]'}    {'ydotdot_{fa} [m/s^2]'}    {'beta'}    {'Tg'}

if useFASTmdl == 1
    % Model Mdl2 from MdlLin folder: 9th order model with state vector:
    % x = [yt zeta xt theta_s yt_dot zeta_dot xt_dot omega_r omega_gr]'
    
    % add path to LinMdl folder and run create linear model for windspeed
    addpath(fullfile(parentDir,'LinMdl'));
    sysOut = createLinearModels(rhoP-3,1);
    
    % Assign matrices
    Av = sysOut.A; Bu = sysOut.B(:,1:2); Cv = sysOut.C(1:4,:);
    Ev = sysOut.B(:,3);
    
    % Initial steady state and output (testX = -inv(Av) * Bv * u0'; testX(1:3)
    x0 = [zeros(1,7),rhoP1(1),rhoP1(1),u1(1),u2(1),d(1)]';
    y0 = x0(8:end-1); %{'omega_r','omega_gr','Tg','beta','V'};
    % State names (input and output names already created)
    sysOut.StateName = {'yt [m]', 'zeta [m]', 'xt [m]', 'theta_s [rad]', ...
        'yt_dot [m/s]', 'zeta_dot [m/s]', 'xt_dot [m/s]',...
        'omega_r [rad/s]', 'omega_gr [rad/s]','Tg [kNm]','beta [rad]'};
    
else
    % System and input matrix
    [AvAll,BvAll] = ABqL_WT_beta(d(1),lambda,beta,Cq);
    
    if idxS == 3
        Av = AvAll(1:idxS,1:idxS);
        Bv = [BvAll(1:idxS,1),AvAll(1:idxS,idxS+1), AvAll(1:idxS,idxS+2)];
    else
        Av = AvAll;
        Bv = BvAll;
    end
    
    Bu = Bv(:,2:3);
    nx = size(Av,2);
    Cv = [zeros(nx-1,1),eye(nx-1)];
    Ev = Bv(:,1);
    
    % Initial steady state and output (testX = -inv(Av) * Bv * u0'; testX(1:3)
    x0 = [0,rhoP1(1),rhoP1(1),u1(1),u2(1),d(1)]';
    y0 = x0(2:end-1);
    sysOut.StateName = {'theta_s [rad]','omega_r [rad/s]',... % State names
        'omega_{gr} [rad/s]','Tg [Nm]','beta [rad]'};
    sysOut.OutputName = [sysOut.StateName(2:end),{'V [m/s]'}]; % Output names
    
    if idxS == 3
        x0 = x0([1:idxS,end]);
        y0 = x0(2:idxS);
        sysOut.StateName = sysOut.StateName(1:idxS);
        sysOut.OutputName = [sysOut.StateName(2:end),{'V [m/s]'}]; % Output names
    end
end

%% Check rank
fprintf('WECS: rank(C*E) is %d; rank(E) is %d\n', rank(Cv*Ev),rank(Ev));
nx = size(Av,1);
ny = size(Cv,1);

% Also create an additional state for the disturbance
tau1 = 1/DT* multTau;
A = [[Av, zeros(nx,1)]; [zeros(1,nx),-tau1]];
C = [Cv, zeros(ny,1)]; %[zeros(1,nx),1]];
E = [Ev;tau1]; %xd_dot = -xd + Ed
B = [Bu; zeros(1,2)];
fprintf('rank(C*E) is %d; rank(E) is %d\n', rank(C*E),rank(E));

%% Calculate observer matrices
H = E*inv((C*E)'*(C*E))*(C*E)'; %#ok<MINV> % inv is fast enough
nx = size(A,1);
T = eye(nx)-H*C;
A1 = T*A;

% Pole placement with Q and R
Qe = [E,B]*[E,B]'; %
Re = C*C';
K1 = lqr(A',C',Qe,Re)'; % lqr not supported in Simulink

% Ham = [A', -C'* Re^(-1) * C; % Hamiltonion matrix -lambda to make it stablizable
%       -Qe, -A];
% [V,D] = eig(Ham); % Divide into stable and unstable eigenvalues
% idxS = diag(D) < 0; % get stable eigenvalues
% Vs = V(:,idxS);% get eigenvectors of stable eigenvalues
% Vus = V(:,~idxS); % not stable
% n11 = sum(idxS);
% P = Vs(n11+1:end,:)/(Vs(1:n11,:));
% K1 = real(Re^(-1) * C * P)';

% Calculate F and K
F = A1 - K1*C; %T*A -K1*C
K = K1 + F*H;
T = eye(size(H*C)) - H*C;

% Sys with no ctrl input: xdot = Ax + Ed, y = Cx
nd = size(E,2); %nu = 0;
nu = size(B,2);
ny = size(C,1);
sys = c2d(ss(A,[E,B],C,zeros(ny,nd+nu)),DT,'zoh'); % ss(A,[E,B],C,zeros(ny,nd+2));

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
ytemp = y0' + cumsum(dy);
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

%% Plot output
figure(1+idxFig)
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


figure(2 + idxFig)
idxEnd = size(x,2)-1;
for idx = 1:idxEnd
    subplot(idxEnd,1,idx)
    plot(t,x(:,idx),'-','Color',cl(1,:));
    hold on;
    plot(t,x1(:,idx),'--','Color',cl(2,:));
    plot(t,xest(:,idx),'k-.');
    axis tight; grid on;
    if idx == 1
        title('WECS states and estimates');
    end
    axis tight;
    grid on;
    ylabel(sysOut.StateName{idx});
end
xlabel('time [s]');

tauMeas = 1/DT;
sysMeasDel = c2d(ss(tf(tauMeas,[1,tauMeas])),DT,'zoh');

%[dXest,~,z] = lsim(sysEst,diffInEnst,t(1:end-1));
windMeas = lsim(sysMeasDel,d,t,d(1)/sysMeasDel.C);

figure(3 + idxFig)
idx = size(x,2);

plot(t, d,'k','Linewidth',1);
hold on;
plot(t,x(:,idx),'-','Color',cl(1,:),'Linewidth',1);
plot(t,xest(:,idx),'--','Color',cl(2,:),'Linewidth',1);
%plot(t,windMeas,'-.','Color',cl(3,:),'Linewidth',2);
legend('wind','wind state','wind est');
ylabel('wind [m/s]')
xlabel('time [s]');
title('Comparison Measured and Estimated Wind');
axis tight;
grid on;
print(gcf,fullfile(figDir,'windEst'), '-dpng');

vafState = 1 - var(x(:,end)- d)/var(d);
vafStateEst = 1 - var(xest(:,end)-d)/var(d);

betaLUT = unique(Lin_points.Pitch); %Select all pitch angle in eq. points
lambdaOP = Lin_points.RSpeed * wecs.rb ./ Lin_points.V; % Calculate TSR in OP
idxLambda = (lambda >= min(lambdaOP) & lambda <= max(lambdaOP)); %index TSR

lambdaLUT = lambda(idxLambda); %select all relevant TSR with index TSR
CqPitch = interp2(lambda,beta,Cq,lambdaLUT,betaLUT);
lambdaInvMat = repmat((lambdaLUT .^ (-2))', size(CqPitch,1),1);
lambdaMat = repmat((lambdaLUT .^ (-2))', size(CqPitch,1),1);
lambdaCqMat = lambdaInvMat .* CqPitch;

lambdaCqVecMin = min(lambdaCqMat(:));
lambdaCqVecMax = max(lambdaCqMat(:));
dLlambdaCqVec = (lambdaCqVecMax - lambdaCqVecMin)/(length(lambdaLUT)-1);
lambdaCqVecLUT = lambdaCqVecMin : dLlambdaCqVec:lambdaCqVecMax;

lambdaTildaMat = lambdaCqVecMax *ones(size(lambdaCqMat));
for idx = 1: length(betaLUT)
    [lambdaCqbeta, idxUnique] = unique(lambdaCqMat(idx,:));
    lambdaTildaMat(idx,:) = interp1(lambdaCqbeta,lambdaLUT(idxUnique),lambdaCqVecLUT,'linear','extrap');
end

%% 3-D Figure
figure(4 + idxFig);

[YLambdaFun,XBeta] = meshgrid(lambdaCqVecLUT,betaLUT);
ax = surf(XBeta,YLambdaFun,lambdaTildaMat);
axis tight;
az = 55;
el = 45;
view(az,el);

xlabel('Pitch \beta [rad]');
ylabel('\lambda^{-2} C_{Q,\beta}(\lambda) [-]');
zlabel('Tip Speed Ratio \lambda [-]')
% title('2D Look-up-Table for \lambda = f(\beta,T_r)');
print(gcf,fullfile(figDir,'nonMonotone3D'), '-dpng');

Jgr = wecs.Jg * wecs.Ng^2;

A = [0,1,-1 ;
    -wecs.Ks/wecs.Jr, -wecs.Bs/wecs.Jr , wecs.Bs/wecs.Jr;
    wecs.Ks/Jgr, wecs.Bs/Jgr, -wecs.Bs/Jgr];
B = [0; 0; -wecs.Ng/Jgr];
E = [0; 1/wecs.Jr;0 ];
C = [0 1 0;0 0 1]; %
fprintf('rank(C*E) is %d; rank(E) is %d\n', rank(C*E),rank(E));

H = E*inv((C*E)'*(C*E))*(C*E)' %#ok<MINV,*NOPTS> display output, inv is fast enough
nx = size(A,1);
T = eye(nx)-H*C;
A1 = T*A;

% Pole placement with Q and R
%Qe = eye(nx); %B*B'; %
%Re = C*C';
Qe = eye(nx); %diag([1,10^11,10^9])*[B,E]*[B,E]'; %eye(nx); %
% r1 = 10^-2;
% r2 = 10^-2;
Re = eye(2); %diag([r1,r2]); %*100;

var1 = var(d)

QXU = diag([10^-4;10^-6;diag([E,B]*[E,B]')]);
QWV = diag([var(d);diag([E,B]*[E,B]');10^(-4);10^(-6)]);

K1 = lqr(A',C',Qe,Re)';

K1a = lqg(ss(A',C', eye(size(A)), zeros(size(C'))),QXU,QWV);

% Calculate F and K
F = A1 - K1*C; %T*A -K1*C
K = K1 + F*H;
T = eye(size(H*C)) - H*C

% Sys with no ctrl input: xdot = Ax + Ed, y = Cx
nd = size(E,2); %nu = 0;
nu = size(B,2);
ny = size(C,1);
sys = ss(A,[E,B],C,zeros(ny,nd+nu)); %c2d(ss(A,[E,B],C,zeros(ny,nd+nu)),DT,'zoh'); % ss(A,[E,B],C,zeros(ny,nd+2));

% Est with with control input: zdot = Fz + TBu + Ky, xhat = z + Hy
sysEstTr = c2d(ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]),DT,'zoh'); %ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]); %

%% Plot lambda for Cqlambda^(-1) for different pitch angles
plotCqlambda2D(Lin_points,wecs,beta,lambda,Cq,cl,idxFig);
print(gcf,fullfile(figDir,'nonMonotone2D'), '-dpng');


%% Plot wind comparison
windSimFile = fullfile('dataOut', 'OutDataWind.mat');
if exist(windSimFile,'file') == 2
    load(windSimFile,'OutDataWind');
else
    return
end

% Get and cut out time vector
timeVec = OutTable.Time;
timeForPlot(1) = 15;
timeForPlot(2) = min(400,timeVec(end));
idxPlot = timeVec >= timeForPlot(1) & timeVec <= timeForPlot(2);
timeVec = timeVec(idxPlot);  % dismiss intiialization

% Set to length of data
idxPlotData = idxPlot(1:length(OutDataWind));
OutDataWindPlot = OutDataWind(idxPlotData,:);
timeVecPlot = timeVec(1:length(OutDataWindPlot));

% wind true for all of relevant times
windTrue = OutDataWindPlot(:,1);
estErr = OutDataWindPlot(:,1) - OutDataWindPlot(:,2);

std(estErr)
1 - var(estErr)/var(windTrue)

%
figure;
subplot(2,1,1)
plot(timeVecPlot, windTrue,'color',cl(1,:));
hold on
plot(timeVecPlot, OutDataWindPlot(:,2),'--','color',cl(2,:));
legend('True','Estimate');
ylabel ('Wind [m/s]')
axis tight; grid on;
subplot(2,1,2)
plot(timeVecPlot,estErr);
legend('True - Estimate');
ylabel ('Wind [m/s]')
xlabel('time [s]')
axis tight; grid on;

print(gcf,fullfile(figDir,'WindEst'), '-dpng');
print(gcf,fullfile(figDir,'WindEst'), '-depsc');

windSimFile = fullfile('dataOut', 'OutDataTr.mat');
if exist(windSimFile,'file') == 2
    load(windSimFile,'OutDataWind');
else
    return
end


estErrTr = OutDataTr(idxPlotData,1) - OutDataTr(idxPlotData,2);
figure;
subplot(2,2,1)
plot(timeVecPlot, OutDataTr(idxPlotData,1),'color',cl(1,:));
hold on
plot(timeVecPlot, OutDataTr(idxPlotData,2),'--','color',cl(2,:));
legend('True','Estimate');
ylabel ('Torque T_r [Nm]')
axis tight; grid on;
subplot(2,2,2)
plot(timeVecPlot,estErrTr);
legend('True - Estimate');
%ylabel ('Wind [m/s]') xlabel('time [s]')
axis tight; grid on;

subplot(2,2,3)
plot(timeVecPlot, OutDataWind(idxPlotData,1),'color',cl(1,:));
hold on
plot(timeVecPlot, OutDataWind(idxPlotData,2),'--','color',cl(2,:));
legend('True','Estimate');
ylabel ('Wind V [m/s]')
xlabel('time [s]')
axis tight; grid on;
subplot(2,2,4)
plot(timeVecPlot,estErr);
legend('True - Estimate');
%ylabel ('Wind [m/s]')
xlabel('time [s]')
axis tight; grid on;

print(gcf,fullfile(figDir,'WindTrEst'), '-dpng');
print(gcf,fullfile(figDir,'WindTrEst'), '-depsc');

% 1 - var(estErrTr)/var(OutDataTr(idxPlotData,1))
