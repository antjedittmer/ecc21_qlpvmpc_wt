%% Prepare workspace
close all; clc %clear;
addpath(genpath(pwd))
parentDir = fileparts(pwd);

%% Define user input
rhoP = 18; %Wind speed: Scheduling parameter, alternative [V,omega,beta]
useFASTmdl = 0; % 0 (Default) 5th order model used, 1 9th order
useFASTwind = 1; % 0 Wind speed step function, 1: NTW Wind from FASTTool
varWind = 1; % Wind variance
multTau = 0.02; %Factor for time constant wind speed low pass ( tau = 1/DT)
multVar = 1; % Change variance NTW

%% Get parameters from Data Dictionary
% The data dictionaries of  WECS test harness with the lookup table (LUT)
LinPointsMFile;
MdlMFile;
LUTMFile;

%% Physical and model constants
DT = 0.008; % sample time
rho = 1.2250;
rb = 63; % m blade radius
Ng = 97; %Gearbox ratio
Jr = (115926 + 3 * 11.776e6) *1; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
Jgr = Jg * Ng^2; % Inertia of the generator in LSS
Ks = 867637000; %Nm/rad Stiffness of the transmission
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
tau_g = 0.01; % for pitch actuator and airspeed

%% Create system matrices A B C E D:
% States: [Theta_s [rad], omega_r [rad/s], omega_gr [rad/s]], u = [Tg]
% d = [Tr(omega_r,beta,V)] y = [omega_r [rad/s], omega_gr [rad/s]]

A = [0,1,-1 ;
    -Ks/Jr, -Bs/Jr , Bs/Jr;
    Ks/Jgr, Bs/Jgr, -Bs/Jgr];
B = [0; 0; -Ng/Jgr];
E = [0; 1/Jr;0 ];
C = [0 1 0;0 0 1]; % 
fprintf('rank(C*E) is %d; rank(E) is %d\n', rank(C*E),rank(E));

H = E*inv((C*E)'*(C*E))*(C*E)' %#ok<MINV,*NOPTS> display output, inv is fast enough
nx = size(A,1);
ny = size(C,1);
T = eye(nx)-H*C
A1 = T*A

% Pole placement with Q and R
Qe = eye(nx); %B*B'; %
Re = C*C';
K1 = lqr(A',C',Qe,Re)';

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
sysEstContTr = ss(F,[T*B,K],eye(nx),[zeros(size(T*B)),H]);

%% To estimate the aerodynamic Torque Tr, use the equation
%      Tr = 0.5 *rho *pi*R^3 *C_Q *V^2
%         = 0.5 *rho *pi*R^3 *C_Q *(omega^2*R^2/lambda^2)
% =>  2* Tr /(rho *pi * R^5 omega^2) = lambda^(-2) C_Q

% Get C_Q,beta as a function of TSR lambda at all relevant pitch angle

[betaLUT,idxLUT] = unique(Lin_points.Pitch); %Select all pitch angle in eq. points
lambdaOP = Lin_points.RSpeed * wecs.rb ./ Lin_points.V; % Calculate TSR in OP
lambdaOPLUT = lambdaOP(idxLUT);
idxLambda1 = 1; %min(lambdaOP) 
idxLambda = (lambda >= idxLambda1 & lambda <= max(lambdaOP)); %index TSR

lambdaLUT = lambda(idxLambda); %select all relevant TSR with index TSR
CqPitch = interp2(lambda,beta,Cq,lambdaLUT,betaLUT);
lambdaInvMat = repmat((lambdaLUT .^ (-2))', size(CqPitch,1),1);
lambdaCqMat = lambdaInvMat .* CqPitch;

lambdaCqFun = nan(length(betaLUT), length(lambdaLUT)); %init lambda^-2 C_Q,beta
CqPitch = nan(length(betaLUT),sum(idxLambda)); % initialize C_Q,beta
for idx = 1: length(betaLUT)
    diffPitch = abs(pitchCt-betaLUT(idx)); % difference to LUT pitch vector
    idxNearestPitch =  min(diffPitch) == diffPitch; % index of nearest pitch
    CqPitch(idx,:) = Cq(idxNearestPitch,idxLambda); % select C_Q,beta
    lambdaCqFun(idx,:) = (lambdaLUT .^ (-2)) .* CqPitch(idx,:)';
end

% Plot tip speed ratio lambda vs. function lambda^-2 C_Q,beta
cl = colormap('lines');
numCl = numel(unique(cl(:,1)));
numPitch = length(betaLUT);
noDiffLineStyles = ceil(numPitch/numCl);
diffLineStyles = {'-','--','-.'};

figure;
for idx = 1: noDiffLineStyles
    nvecAll = (1:numCl) + (idx -1)*numCl;
    nvec = nvecAll(nvecAll<= numPitch);
    plot(lambdaLUT,lambdaCqFun(nvec,:),'LineStyle',diffLineStyles{idx});
    hold on;
end

X = [1, min(lambdaOP), min(lambdaOP),1];
Y = [0,0, max(lambdaCqMat(:)),  max(lambdaCqMat(:))];
fill(X,Y,'k','FaceAlpha',0.2,'LineStyle','none');
axis tight; grid on

xlabel('Tip Speed Ratio \lambda [-]');
ylabel('\lambda^{-2} C_{Q,\beta}(\lambda)[-]');

% 
% idxSpacing = 1;
% for idx = 1: idxSpacing: length(lambdaOPLUT)
%     idxVec = abs(lambdaLUT - lambdaOPLUT(idx));
%     idxOP = idxVec  == min(idxVec)
%     plot(lambdaOPLUT(idx),lambdaCqFun(idx,idxOP),'o','MarkerFaceColor',cl(idx,:),'MarkerEdgeColor',cl(idx,:) );
% end


cellLeg = regexp(sprintf('Pitch: %2.2f [deg],', betaLUT(1:idxSpacing:end)*180/pi),',','split');
legend(cellLeg(1:end-1));
%title('Monotonically non-increasing function \lambda^{-2}C_{Q,\beta}(\lambda)');

figure; plot(Lin_points.V, lambdaOP)

axis tight; grid on
xlabel('Wind speed [m/s]');
ylabel('\lambda[-]');


%% Check if wind speed can be deduced from lambda^(-2) C_Q
vLut = unique(Lin_points.V);
lambdaCqPlot = nan(length(vLut),1);
for idx = 1: length(vLut)
    lambdaTemp = lambdaOP(idx);
    betaTemp = Lin_points.Pitch(idx);
    CQtemp = interp2(lambda,pitchCt,Cq,lambdaTemp,betaTemp);
    lambdaCqPlot(idx) = lambdaTemp^(-2) *CQtemp;
end

figure; plot(vLut, lambdaCqPlot)

%%
%% Set-up wind test signal
% Get time and input vectors to simulate with lsim

t = 0:DT:100; % time vector

dataInPath = fullfile(parentDir ,'dataIn');
load(fullfile(dataInPath,'OutDataWind18NTW.mat'),'OutTable');
%load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');

multVar = 0;
varWind = 0;

windTest = OutTable.Wind1VelX;
meanWindTest = mean(windTest); % varWindTest = var(windTest);
OutTable.Wind1VelX = multVar*(windTest - meanWindTest) + meanWindTest;

% t = OutTable.Time;
V = OutTable.Wind1VelX(1:length(t));

d1 = varWind*randn(size(V));
d = 6*ones(size(V)); %  + d1;  % wind input load('d'); d = d';


%% Set up control inputs as feedforward  and initial state rotor speed
% u = [Tg_ref beta_ref V]'

% Get torque, pitch and and rotor speed corresponding to wind input
Tr = nan(size(d));
Tg = nan(size(d));
for idx = 1: length(d)
    idxVec = abs(Lin_points.V - d(idx)); % difference to wind speed
    idxNearestV =  min(idxVec) == idxVec; % index of nearest wind speed
    Tg(idx) = Lin_points.Torque(idxNearestV )*10^3; % GenTorque [Nm]
    % lambda, beta and Cq at current wind speed
    lambdaTemp = lambdaOP(idxNearestV);% TSR lambda [-]
    betaTemp = Lin_points.Pitch(idxNearestV);% Pitch [rad]
    CQtemp = interp2(lambda,pitchCt,Cq,lambdaTemp,betaTemp); %Cq(lambda, pitch)
    
    Tr(idx) = 0.5* rho *pi * rb^3 * CQtemp *d(idx)^2; %Aerodynamic torque
        
end

idx1 = abs(Lin_points.V - d(1)); % difference to wind speed
idxNearestV1 =  min(idx1) == idx1;
x0 = [0,Lin_points.RSpeed(idxNearestV1),Lin_points.RSpeed(idxNearestV1)]';

in1 = [Tr,Tg];
diffIn = diff(in1);

[yhat,~,xhat] = lsim(sys,in1,t,x0);
y1 = yhat;% + y0';
x1 = xhat;%+ x0';


inEst = [Tg,y1];
%
z0 = T*x0; %(eye(size(H*C))- H*C);
[xesthat,~,zhat] = lsim(sysEstTr,inEst,t,z0);
z1 = zhat;% + z0';
xest1 = xesthat;% + x0';

  
% State and output names
sysOut.StateName = {'theta_s [rad]','omega_r [rad/s]',...
    'omega_{gr} [rad/s]'};
sysOut.OutputName = [sysOut.StateName(2:3)];

%
figure
cl = colormap('lines');

ny = size(y1,2);
for idx = 1:ny
    subplot(ny, 1,idx)
    plot(t,y1(:,idx), 'Color',cl(1,:),'Linewidth',1);
    if idx == 1
        title('WECS outputs');
    end
    axis tight;
    grid on;
    ylabel(sysOut.OutputName{idx});
end
xlabel('time [s]');

figure
idxEnd = size(x1,2);
for idx = 1:idxEnd
    subplot(idxEnd,1,idx)
    plot(t,x1(:,idx),'-','Color',cl(1,:));
    hold on;
    plot(t,xest1(:,idx),'k--');
    axis tight; grid on;
    if idx == 1
        title('WECS states and estimates');
    end
    axis tight;
    grid on;
    ylabel(sysOut.StateName{idx});
end
xlabel('time [s]');

% J*omega_r_hat_dot = -Ks*Theta_s - Bs*omega_r_hat + Bs*omega_gr_hat + Tr
% => Tr = Jr*omega_r_hat_dot + Ks*Theta_s + Bs*omega_r_hat - Bs*omega_gr_hat
sysEstAug = ss(sysEstTr.A,sysEstTr.B,[sysEstTr.C;sysEstTr.A(2,:)],[sysEstTr.D;sysEstTr.B(2,:)]);
[xesthatAug,~,zhatAug] = lsim(sysEstAug,inEst,t,z0);

omegaRhat_dot = xesthatAug(:,4);
omegaRhat = xesthatAug(:,2);
omegaGhat = xesthatAug(:,3);
ThetaShat = xesthatAug(:,1);

TgNg = Ng * Tg;


% A = [0,1,-1 ;
%     -Ks/Jr, -Bs/Jr , Bs/Jr;
%     Ks/Jgr, Bs/Jgr, -Bs/Jgr];
% 
% B = [0; 0; -Ng/Jgr];
% E = [0; 1/Jr;0 ];
% C = [0 1 0;0 0 1];

Trhead = Jr*omegaRhat_dot + Bs*omegaRhat + ....
    -Bs*omegaGhat  + Ks*ThetaShat;
%TgNg;
    
    %  ;

figure;
subplot(2,1,1)
plot(t,Tr,t,Trhead);
axis tight; grid on;
subplot(2,1,2)
plot(t,Tr-Trhead);
xlabel('time'); 



%Jr = inv(sys.B(2,1));
sys.B(2,2)

% Trest = nan(length(t),1);
% for idx = 2:1:length(t)
% Trest(idx-1) = inv(sys.B(2,1))*(xest1(idx-1,2)-sys.A(2,:)*xest1(idx-1,:)'-sys.B(2,2)*Tg(idx-1));
% end
% figure; plot(t,Tr,t,Trest)

% 2 Qa ? ? R5 ?r
%% Calculate lambdaCqFun
lambdaCqFun = 2 * Trhead ./(rho * pi *rb^5 * xest1(:,2).^2);

% without known beta -> directly V
Vest1 = interp1(lambdaCqPlot,vLut,lambdaCqFun);
figure; plot(t,d', t,Vest1');




