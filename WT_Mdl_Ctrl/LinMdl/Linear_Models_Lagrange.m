clear; clc; close all
% addpath('C:\Program Files\export_fig'); %C:\Users\bs\Desktop\export_fig-master\export_fig-master;
% addpath C:\Werk_Data\Software\hline_vline;
% load 'NREL5MW_linearised_reduced'
load 'NREL5MW_linearised'
load 'NREL5MW_CPdata_new3'
Lin_points = Lin;

%% Plot parameters
% yt is the axial displacement of the nacelle, $\zeta$ is the angular blade displacement out of the plane of rotation and ?r and ?g are the angular positions of the
% X=[yt zeta thetha_r thetha_gr yt_dot zeta_dot omega_r omega_gr ]'
% uact=[Tg_ref beta_ref V ]'
% y=[ yt_dot zeta_dot omega_r omega_gr]'
% u=[beta_ref wg_ref]'
% d=v
% y=[yt_dot zeta_dot omega_r omega_z]'
createBodePlots = 1;
plotVisible = 'on';
% sysOutputname=[{'omega [rad/s]'}; {'\ddot{x}_{fa} [m/s^2]'}; {'\ddot{x}_{ss} [m/s^2]'}];
sysOutputname =[{'xdotdot_{fa} [m/s^2]'};{'rb*zeta[m]'}; {'omega_r [rad/s]'};];
sysInputname =[{'T_g [kNm]'}; {'\beta [rad]'};{'V [m/s]'} ];

%% Create output folder
figFolder = 'FigureLagrange_201014'; %'Figures';
if ~isfolder(figFolder)
    mkdir(figFolder);
end

%% Define parameters
% 5MW FAST WIND TURBINE (NREL/TP-500-38060)Definition of a 5-MW Reference Wind Turbine for Offshore
%System Development
% Cut-In, Rated, Cut-Out Wind Speed 3 m/s, 11.4 m/s, 25 m/s
% Cut-In, Rated Rotor Speed 6.9 rpm, 12.1 rpm
%wr,wg in rad/s, Tg in kN, Pg in kW,
Ng = 97; %Gearbox ratio
mt = 347460 + 240000;%347460;%kg Mass of the tower and nacelle
mb = 17740; %kg Mass of each blade
Jr = 115926 + 3 * 11.776e6; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator
wnt = 2.1; %2.1 ; %2.9003; %Natural frequency of tower
wnb = 2.0205; %Natural frequency of blade
zetat = 0.01; %damping ratio of tower
zetab = 0.00477465;%damping ratio of blade

% K=wn^2
Kt = wnt^2 * mt; % Stiffness of the tower s^2 + B/Ms + K/m
Kb = wnb^2 * mb; %Stiffness of each blade

Ks = 867637000; %Nm/rad Stiffness of the transmission
% 2*zeta*wn=B
Bt = 2*zetat*wnt *mt; % Damping of the tower
Bb = 2*zetab*wnb *mb ; %Damping of the blade
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission

N = 3;%Number of blades
rb = 63; % m blade radius
rho = 1.225; % Air density (kg/m^3)
Bg = 0.9; % Tg= Bg(wg - wz).
tau = 0.1; % for pitch actuator
kappa = 0.01;

%% Parameters for tower model
% kfa = 1; %should this be a plus or a min?
% Mfa = 4.38e5; % mt+N*mb: 400680
% Cfa = 2*6421;
% Kfa = 1.846e6;
% 
% % mb =(Mfa - mt)/N;
% mt = Mfa;
% Kt = Kfa;
% Bt = Cfa; % 2*zetat*wnt *mt

% Kb = wnb^2*mb; %
% Bb  = 2*zetab*wnb*mb;

%% Model
M =[mt+N*mb N*mb*rb 0 0;
    N*mb*rb N*mb*rb^2 0 0;
    0 0 Jr 0;
    0 0 0 Jg*Ng^2];

Ce = [Bt 0 0 0;
    0 N*Bb*rb^2 0 0;
    0 0 Bs -Bs;
    0 0 -Bs Bs];

K = [Kt 0 0;
    0 N*Kb*rb^2 0;
    0 0 Ks;
    0 0 -Ks];

Q = [1 0 0 ;
    rb 0 0 ;
     0 1 0;
     0 0 -Ng];
L = [eye(3), [0;0;-1]];
A1 = [zeros(3), L;...
    -M\K, -M\Ce];
B1=[zeros(3); M\Q];

vecLambda = 1:150;
Rotor_Lambda = Rotor_Lamda(vecLambda);
Rotor_Pitch = Rotor_Pitch * pi/180; % Rotor Pitch deg ->rad
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');

%% Compute gradients

[Cqdlambda,Cqdbeta] = gradient(Cq,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Ctdlambda,Ctdbeta] = gradient(Ct,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
Clambda = Cq./Xq.^2;
[Clambdadlambda,Clambdadpitch] = gradient(Clambda,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));

%% Initialize Wnn ZZ PP matrices
sizeA = 5; % number of states, ToDo AD: Is there an automatic way to determine this?
Wnn = NaN(sizeA,length(Lin.V));
ZZ = NaN(sizeA,length(Lin.V));
PP = NaN(sizeA,length(Lin.V));

%% Options for Bode plots
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz'; opts.MagUnits ='Abs'; opts.MagScale ='log'; opts.XLim ={[1e-2 1e1]};

%% Initialize name of models and load
idxModelNames = 1:2; %1:5
modelNames =[{'NREL5MW_linearised'}; {'NREL5MW_linearised_Gen_Fa_Flap1'};
    {'NREL5MW_linearised_Gen_Fa_SS_Flap1'};
    {'NREL5MW_linearised_Gen_DriveTrain_Fa_SS_Flap1'}; {'NREL5MW_linearised'} ];
modelNames = modelNames(idxModelNames);

modelsL = [{'FAST 30 states'};
    {'FAST 10 states'};
    {'Fast Gen+FA+SS+1^{st} flap mode (12 states)'};
    {'Fast Gen+drive train+FA+SS+1^{st} flap mode (14 states)'};
    {'Fast Gen+FA (4 states)'}];

modelsL = modelsL(idxModelNames);

lenModelNames = length(modelNames);
sysmCell = cell(lenModelNames,1);
sysmCell_MBC = cell(lenModelNames,1);
for idxModels = 1: lenModelNames
    load(modelNames{idxModels},'sysm','sysm_MBC');
    sysmCell{idxModels} = sysm;
    sysmCell_MBC{idxModels} = sysm_MBC;
end

%% Initialize natural frequencies, damping, pole vectros
lenLinV = length(Lin.V);
sizeA = size(sysm{1}.A,1);
Wn = cell(lenLinV,sizeA);
Z = cell(lenLinV,sizeA);
P = cell(lenLinV,sizeA);
sizeA_MBC = size(sysm_MBC{1}.A,1);
Wnmbc = cell(lenLinV,sizeA);
Zmbc = cell(lenLinV,sizeA);
Pmbc = cell(lenLinV,sizeA);


%% Check original Output names 
sysOutNamesOrig = sysm{1}.OutputName([12 6 38]);
disp(['FAST sysOutputnames: ', sprintf('%s, ', sysOutNamesOrig{:})])
disp(['TU Delft sysOutputnames: ', sprintf('%s, ', sysOutputname{:})])

sysInNamesOrig = sysm{1}.InputName([8 9 1]);
disp(['FAST sysOutputnames: ', sprintf('%s, ', sysInNamesOrig{:})])
disp(['TU Delft sysOutputnames: ', sprintf('%s, ', sysInputname{:})])


%% Loop over all air speeds in LinV
for index = 1 %[1,length(Lin.V)]%1 : length(Lin.V) %[1,length(Lin.V)] % [1,2,length(Lin.V)-1,length(Lin.V)]
    
    %% Variables dependent on wind speed V
    Vbar = Lin_points.V(index);
    Lambda_bar=Lin_points.RSpeed(index)*rb/Vbar;
    omegabar=Lin_points.RSpeed(index);
    beta_bar=Lin_points.Pitch(index); %/pi*180;
    
    dlambdadomega = rb/Vbar;
    dlambdadV = -omegabar*rb/Vbar^2;
    dlambdadxfa = omegabar*rb/Vbar^2;
    
    %% Grid points
    dCqdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
    dCtdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdlambda,Lambda_bar,beta_bar,'cubic');
    dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
    %dCqdbetab = dCqdbetabDeg * 180/pi;
    dCtdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic');
    Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
    Ctb = interp2(Rotor_Lambda,Rotor_Pitch,Ct,Lambda_bar,beta_bar,'cubic');
    
    %% Ft=0.5*pi*rb^2*Ct*v^2
    % bt = -0.5*rho*pi*rb^2*dCtdlambdab*dlambdadomega*Vbar^2; : AD should
    % this be +
    % Ktv = 0.5*rho*pi*rb^2*(dCtdlambdab*dlambdadV*Vbar^2+2*Vbar*Ctb);%dFt/dV
    % Ktb = 0.5*rho*pi*rb^2*dCtdbetab*Vbar^2;%dFt/dbeta
    dFtdomega = 0.5*rho*pi*rb^2*dCtdlambdab*dlambdadomega*Vbar^2;%dFt/domega
    dFtdV = 0.5*rho*pi*rb^2*(dCtdlambdab*dlambdadV*Vbar^2+2*Vbar*Ctb);%dFt/dV
    dFtdbeta = 0.5*rho*pi*rb^2*dCtdbetab*Vbar^2;%dFt/dbeta
    
    %% Tr=0.5*pi*rb^3*Cq*v^2
    %br = -0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;%-dTr/domega
    %Krv = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadV*Vbar^2+2*Vbar*Cqb);%dTr/dV
    %Krb = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;%dTr/dbeta
    dTrdomega = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;%dTr/domega
    dTrdV = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadV*Vbar^2+2*Vbar*Cqb);%dTr/dV
    dTrdbeta = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;%dTr/dbeta
    
  
    %% Linear Model
    % X=[yt zeta thetha_r thetha_gr yt_dot zeta_dot omega_r omega_gr beta ]'
    % For stability
    % X=[yt zeta thetha_s yt_dot zeta_dot omega_r omega_gr beta]'
    % y=[ yt_ddot rb*zeta_dot omega_r]'
    % uact=[Tg_ref beta_ref V ]'
    % Ve = V - yt_dot [m/s] - rb [m]  * zeta_dot [1/s]
    
    A = [A1(1:3,:), zeros(3,2);...
        A1(4,1:3),A1(4,4)-B1(4,1)*dFtdV,A1(4,5)-B1(4,1)*dFtdV*rb, A1(4,6)+B1(4,1)*dFtdomega ,A1(4,7),0, B1(4,1)*dFtdbeta;...
        A1(5,1:3),A1(5,4)-B1(5,1)*dFtdV,A1(5,5)-B1(5,1)*dFtdV*rb, A1(5,6)+B1(5,1)*dFtdomega ,A1(5,7),0, B1(5,1)*dFtdbeta;...
        A1(6,1:3),A1(6,4)-B1(6,2)*dTrdV,A1(6,5)-B1(6,2)*dTrdV*rb, A1(6,6)+B1(6,2)*dTrdomega ,A1(6,7),0, B1(6,2)*dTrdbeta;...
        A1(7,1:7), B1(7,3), 0;...
        zeros(1,7),  -1/kappa, 0;...
        zeros(1,8),  -1/tau];
    
    B = [zeros(3);...
       zeros(1,2),B1(4,1)*dFtdV;...
       zeros(1,2),B1(5,1)*dFtdV;...
       zeros(1,2),B1(6,2)*dTrdV;...
       zeros(1,3);...
       1/kappa, 0,0;...
       0,1/tau,0];
    C = [A(4,:);...
        0     rb     zeros(1,7);...
        zeros(1,5),     1     zeros(1,3)];
    D = [B(4,:);zeros(2,3)];
        
    %% Create a state-space model
    sys = ss(A,B*diag([1000,1,1]),C,D); % to give input Tg in kN/m
    sys.Outputname = sysOutputname;
    sys.Inputname = sysInputname;
    
    sysCell{index} = sys;
    Torque_g_bar = 0.5*rho*pi*rb^3*Cqb*Vbar^2/Ng; % Nm
    u_bar_vec(:,index) = [Torque_g_bar, beta_bar * pi/180,Vbar]';
    x_bar_vec(:,index)  = - sys.A\sys.B *u_bar_vec(:,index); % X(s) = (sI - A)^1 *B *U(s) (xdot = 0 = Ax +Bu)
    y_bar_vec(:,index)  = sys.C * x_bar_vec(:,index) + sys.D * u_bar_vec(:,index);

    
    %% Make a Bode plot
    if createBodePlots == 1
        ff = figure(index*100+3);%ff.Units='centimeters';
        ff.Position = [520   267   560   531];
        ff.Color = 'white'; ff.Visible = plotVisible;
        
        bodemag(sys,'k'); hold on; grid on;
            
        for indexModels=1:lenModelNames
            sysm =  sysmCell{indexModels};
            model = [1 0 0; 0 1 0; 0 0 2*pi/60]*sysm{index}([12 6 38],[8 9 1])*[1000 0 0; 0 1 0; 0 0 1]; %pi/180
            model.OutputName = sysOutputname;
            model.InputName = sysInputname;
            bodemag(model,opts); hold on
            grid on
        end
        
        legCell = [{'Bianchi 8 states'};modelsL];
        cl = colormap('lines');
        strExtr = '';
        for ilegC = 1: length(modelsL)
            strExtr = [strExtr,...
                ', \color[rgb]{',num2str(cl(ilegC,:)),'}',modelsL{ilegC}]; %#ok<AGROW> This fast
        end
        
        %legend(legCell,'location','southeast')
        
        titleStr = sprintf('Bode mag. plot for V = %d, pitch = %2.2f, and TSR = %2.4f',Vbar,beta_bar,Lambda_bar);
        titleCell = {titleStr, [legCell{1},strExtr]};
        title(titleCell)
        
        figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),'.png'];
        figFolderStr = fullfile(figFolder,figStr);
        % export_fig(figFolderStr,'-dpng','-m2','-nocrop') %ToDo AD Check why
        print(figFolderStr, '-dpng');
    end
    
end



return; %
ff = figure(50);ff.Units='centimeters'; ff.Position=[2 2 45 27]; ff.Color='white';
subplot(2,3,4)
Damping=[Zmbc{:,4}];
Freq=[Wnmbc{:,4}]/2/pi;
plot(Lin.V,Damping(3:end,:)); hold on
ylabel('Modal damping')
xlabel('Wind Speed [m/s]')
subplot(2,3,1)
plot(Lin.V,Freq(3:end,:)); hold on
ylabel('Frequency in Hz')
xlabel('Wind Speed [m/s]')
title('Modal frequencies and damping vs V [MBC transformed]')
ylim([0 3])

subplot(2,3,5)
Damping = [Z{:,4}] ;
Freq= [Wn{:,4}]/2/pi;
plot(Lin.V,Damping(3:end,:)); hold on
ylabel('Modal damping')
xlabel('Wind Speed [m/s]')
subplot(2,3,2)
plot(Lin.V,Freq(3:end,:)); hold on
ylabel('Frequency in Hz')
xlabel('Wind Speed [m/s]')
title('Modal frequencies and damping vs V ')
ylim([0 3])
subplot(2,3,6)
Damping = ZZ;
Freq = Wnn/2/pi;
plot(Lin.V,Damping(3:end,:)); hold on
ylabel('Modal damping')
xlabel('Wind Speed [m/s]')
subplot(2,3,3)
plot(Lin.V,Freq(3:end,:)); hold on
ylabel('Frequency in Hz')
xlabel('Wind Speed [m/s]')
title('Modal frequencies and damping vs V ')
ylim([0 3])

figStr = 'Figure_Damping_Freq.png';
figFolderStr = fullfile(figFolder,figStr);
%export_fig(figStr,'-dpng','-m2','-nocrop')
print(figFolderStr, '-dpng');




