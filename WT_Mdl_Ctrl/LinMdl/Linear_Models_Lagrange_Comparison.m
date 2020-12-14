function wecs = Linear_Models_Lagrange_Comparison(speedVec,useTUDelftParam,figFolder, useActuatorStates,figNoAdd)

%% Set defaults
if nargin == 0
    speedVec = [1,10,22];
    useTUDelftParam = 0;
    figFolder = 'FigureLagrange_201015_comp';
    useActuatorStates = 0;
    figNoAdd = 0;
    close all;
end

%% Initialize: Set path
% addpath('C:\Program Files\export_fig'); %C:\Users\bs\Desktop\export_fig-master\export_fig-master;
% addpath C:\Werk_Data\Software\hline_vline;
% load 'NREL5MW_linearised_reduced' % Not available
addpath(fullfile(fileparts(pwd),'dataIn'));
%load('NREL5MW_linearised','sysm', 'Lin', 'sysm_MBC')
% load('NREL5MW_linearised_4to25','sysm', 'Lin')
load('NREL5MW_CPdata_new3','Rotor_Lamda', 'Rotor_Pitch', ...
    'Rotor_cQ', 'Rotor_cT');

%% Initialize name of models and load into workspace
idxModelNames = 1:3 %1:2; %1:5
modelNames =[ {'NREL5MW_linearised'}; %{'NREL5MW_linearised_4to25'};
    {'NREL5MW_linearised_Gen_Fa_Flap1'};
    {'NREL5MW_linearised_Gen_Fa_SS_Flap1'};
    {'NREL5MW_linearised_Gen_DriveTrain_Fa_SS_Flap1'}; {'NREL5MW_linearised'} ];
modelNames = modelNames(idxModelNames);

modelsL = [{'FAST 30 states'};
    {'FAST 10 states'};
    {'Fast Rot+FA+SW+1^{st} flap mode (12 states)'};
    {'Fast Rot+drive train+FA+SW+1^{st} flap mode (14 states)'};
    {'Fast Rot+FA (4 states)'}];

modelsL = modelsL(idxModelNames);

lenModelNames = length(modelNames);
sysmCell = cell(lenModelNames,1);
sysmCell_MBC = cell(lenModelNames,1);

idxModels = 1;

load(modelNames{idxModels},'sysm','Lin'); %'sysm_MBC');
sysmCell{idxModels} = sysm;
for idxModels = 2 : lenModelNames
    load(modelNames{idxModels},'sysm','sysm_MBC');
    sysmCell{idxModels} = sysm;
    sysmCell_MBC{idxModels} = sysm_MBC;
end


%% Plot parameters
% yt is the axial displacement of the nacelle, $\zeta$ is the angular blade displacement
% out of the plane of rotation and $\omega_r$ and $\omega_r$ are the angular positions of the
% X=[yt zeta thetha_r thetha_gr yt_dot zeta_dot omega_r omega_gr ]'
%uact=[Tg_ref beta_ref V ]'
%y=[ yt_dot zeta_dot omega_r omega_gr]'
%u=[beta_ref wg_ref]'
% d=v
% y=[yt_dot zeta_dot omega_r omega_z]'
createBodePlots = 1;
plotVisible = 'on';
% sysOutputname=[{'omega [rad/s]'}; {'\ddot{x}_{fa} [m/s^2]'}; {'\ddot{x}_{ss} [m/s^2]'}];
% sysOutputname=[{'xdotdot_{fa} [m/s^2]'};{'rb*zeta[m]'}; {'omega_r [rad/s]'};];
sysOutputname = [{'omega_r [rad/s]'};{'ydotdot_{fa} [m/s^2]'}];
sysInputname =[{'Torque [kNm]'}; {'\beta [rad]'};{'V [m/s]'}];

%% Create output folder
if ~isfolder(figFolder)
    mkdir(figFolder);
end

%% Define parameters
% 5MW FAST WIND TURBINE (NREL/TP-500-38060) Definition of a 5-MW Reference 
% Wind Turbine for Offshore System Development
% Cut-In, Rated, Cut-Out Wind Speed 3 m/s, 11.4 m/s, 25 m/s
% Cut-In, Rated Rotor Speed 6.9 rpm, 12.1 rpm
% wr,wg in rad/s, Tg in kN, Pg in kW,
wecs.Ng = 97; %Gearbox ratio
wecs.mh =  56780; % hub mass
wecs.mb = 17740; %kg Mass of each blade
wecs.mn = 240000; % Nacellle mass
wecs.mtower = 347460;%347460;%kg Mass of the tower and nacelle
wecs.mt = 0.25*wecs.mtower+wecs.mn+wecs.mh+3*wecs.mb;
wecs.H =  87.6;

wecs.Jr = 115926 + 3 * 11.776e6; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
wecs.Jrg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator
f0 = 0.324;
f0sw = 0.3120;
wecs.wnb = 0.6993 * 2*pi; %2.0205; %Natural frequency of blade
wecs.wnt = f0 * 2*pi; % wecs.wnb; 0.3240
wecs.zetat = 1/100; %damping ratio of tower (Table 6.2)
wecs.zetab = 0.477465;% damping ratio of blade /100
% K=wn^2

%if ~useTUDelftParam
wecs.Kt = wecs.wnt^2 * wecs.mt; % Stiffness of the tower s^2 + B/Ms + K/m
wecs.Bt = 2 *wecs.zetat * wecs.wnt *wecs.mt; % tower 2*6421;
%end
wecs.Kb = wecs.wnb^2 * wecs.mb; %Stiffness of each blade
wecs.Bb = 2 *wecs.zetab*wecs.wnb *wecs.mb ; %Damping of the blade
wecs.Ks = 867637000; %Nm/rad Stiffness of the transmission
% 2*zeta*wn=B
wecs.Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission

wecs.N = 3;%Number of blades
wecs.rb = 63; % m blade radius
rho = 1.225; % Air density (kg/m^3)
% Bg = 0.9; % Tg= Bg(wg - wz). Unused because we use Tg as input
tau = 0.1; % for pitch actuator
kappa = 0.01;

% Tower equivalent mass, MT 438,000 kg
% Tower equivalent damping, CT 6421 wecs.N s/m
% Tower equivalent stiffness, KT 1,846,000 wecs.N/m

%% Parameters for tower model
% wecs.mt = 4.38e5; % mt+wecs.N*mb: 400680
% wecs.Bt = 12842; %2*6421;
% wecs.Kt = 1.846e6;
wecs.Kt = wecs.Kt;
wecs.Bt = wecs.Bt; %wecs.Bt;
Kss = wecs.Kt;
% wecs.wnt = sqrt(wecs.Kt/wecs.mt); 

%% Lagrange's Model matrices
M =[wecs.mt wecs.N*wecs.mb*wecs.rb 0 0; % M(1,1) mt+wecs.N*mb  -> mt (already includes blades)
    wecs.N*wecs.mb*wecs.rb wecs.N*wecs.mb*wecs.rb^2 0 0;
    0 0 wecs.Jr 0;
    0 0 0 wecs.Jrg*wecs.Ng^2];

Ce = [wecs.Bt 0 0 0;
    0 wecs.N*wecs.Bb*wecs.rb^2 0 0;
    0 0 wecs.Bs -wecs.Bs;
    0 0 -wecs.Bs wecs.Bs];

K = [wecs.Kt 0 0;
    0 wecs.N*wecs.Kb*wecs.rb^2 0;
    0 0 wecs.Ks;
    0 0 -wecs.Ks];

Q = [1 0 0 ;
    wecs.rb 0 0 ;
    0 1 0;
    0 0 -wecs.Ng];

L = [eye(3), [0;0;-1]];

A1 = [zeros(3), L;...
    -M\K, -M\Ce];

B1 = [zeros(3); M\Q];

%% Select relevant section LUT
vecLambda = 1:150;
Rotor_Lambda = Rotor_Lamda(vecLambda);
Rotor_cQ = Rotor_cQ(:,vecLambda);
Rotor_cT = Rotor_cT(:,vecLambda);
Rotor_Lamda = Rotor_Lambda;
Rotor_Pitch = Rotor_Pitch * pi/180; % Rotor Pitch deg ->rad
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');

%% Compute gradients LUT
[Cqdlambda,Cqdbeta] = gradient(Cq,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
[Ctdlambda,Ctdbeta] = gradient(Ct,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));

%% Initialize Wnn ZZ PP matrices
sizeA = 5; % number of states, ToDo AD: Is there an automatic way to determine this?
Wnn = NaN(sizeA,length(Lin.V));  %#ok<*NASGU>
ZZ = NaN(sizeA,length(Lin.V));
PP = NaN(sizeA,length(Lin.V));

%% Options for Bode plots
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz'; opts.MagUnits ='Abs'; opts.MagScale ='log'; opts.XLim ={[1e-2 1e1]};

%% Initialize natural frequencies, damping, pole vectros
lenLinV = length(Lin.V);
sizeA = size(sysm{1}.A,1);

P = cell(lenLinV,sizeA);
if exist('sysm_MBC','var') == 1
    sizeA_MBC = size(sysm_MBC{1}.A,1);
end
Pmbc = cell(lenLinV,sizeA);
sysCell = cell(lenLinV,1);
u_bar_cell = cell(lenLinV,1);
x_bar_cell =  cell(lenLinV,1); % X(s) = (sI - A)^1 *B *U(s) (xdot = 0 = Ax +Bu)
y_bar_cell =  cell(lenLinV,1);

% Wn = cell(lenLinV,sizeA);
% Z = cell(lenLinV,sizeA);
% Wnmbc = cell(lenLinV,sizeA);
% Zmbc = cell(lenLinV,sizeA);


%% Check original Input and Output names
idxFASTOutput = [38 12] ; % 38 12 [12 6 38]
idxFASTInput = [8 9 1];

sysOutNamesOrig = sysm{1}.OutputName([idxFASTOutput,6]);
disp(['TU Delft sysOutputnames: ', sprintf('%s, ', sysOutputname{:})])
disp('FAST sysOutputnames:')
fprintf('%s\n', sysOutNamesOrig{:});
fprintf('\n')
disp(['TU Delft sysInputnames: ', sprintf('%s, ', sysInputname{:})])
sysInNamesOrig = sysm{1}.InputName(idxFASTInput);
disp('FAST sysInputnames: ');
fprintf('%s\n', sysInNamesOrig{:})


%% Loop over all air speeds in LinV
for index =  speedVec%[1,9,10,length(Lin.V)]%1 : length(Lin.V) %[1,length(Lin.V)] % [1,2,length(Lin.V)-1,length(Lin.V)]
    
    %% Variables dependent on wind speed V
    Vbar = Lin.V(index);
    Lambda_bar = Lin.RSpeed(index)*wecs.rb/Vbar;
    omegabar = Lin.RSpeed(index);
    beta_bar = Lin.Pitch(index); %/pi*180;
    
    dlambdadomega = wecs.rb/Vbar;
    dlambdadV = -omegabar*wecs.rb/Vbar^2;
    dlambdadxfa = omegabar*wecs.rb/Vbar^2; %s/m
    dlambdadxdotfa = -1*omegabar*wecs.rb/Vbar^2; %s/m   xdotfa m/s
    
    %% Grid points
    dCqdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
    dCtdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdlambda,Lambda_bar,beta_bar,'cubic');
    dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
    dCtdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic'); %dCqdbetab = dCqdbetabDeg * 180/pi;
    Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
    Ctb = interp2(Rotor_Lambda,Rotor_Pitch,Ct,Lambda_bar,beta_bar,'cubic');
    
    %% Ft=0.5*pi*rb^2*Ct*v^2
    dFtdomega = 0.5*rho*pi*wecs.rb^2*dCtdlambdab*dlambdadomega*Vbar^2;%dFt/domega
    dFtdV = 0.5*rho*pi*wecs.rb^2*(dCtdlambdab*dlambdadV*Vbar^2+2*Vbar*Ctb);%dFt/dV
    dFtdbeta = 0.5*rho*pi*wecs.rb^2*dCtdbetab*Vbar^2;%dFt/dbeta
    
    %% Tr=0.5*pi*rb^3*Cq*v^2
    dTrdomega = 0.5*rho*pi*wecs.rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;%dTr/domega
    dTrdV = 0.5*rho*pi*wecs.rb^3*(dCqdlambdab*dlambdadV*Vbar^2+2*Vbar*Cqb);%dTr/dV
    dTrdbeta = 0.5*rho*pi*wecs.rb^3*dCqdbetab*Vbar^2;%dTr/dbeta
    
    %% Linear Model
    % X =[yt zeta Thetha_s yt_dot zeta_dot omega_r omega_gr beta]'
    % y=[ yt_ddot rb*zeta_dot omega_r]' -> y=[ yt_ddot rb*zeta_dot omega_r]'
    % uact=[Tg_ref beta_ref V]'
    % Ve = V - yt_dot [m/s] - rb [m]  * zeta_dot [1/s]
    
    A9DoF = [A1(1:3,:), zeros(3,2);...
        A1(4,1:3),A1(4,4)-B1(4,1)*dFtdV,A1(4,5)-B1(4,1)*dFtdV*wecs.rb, A1(4,6)+B1(4,1)*dFtdomega ,A1(4,7),0, B1(4,1)*dFtdbeta;...
        A1(5,1:3),A1(5,4)-B1(5,1)*dFtdV ,A1(5,5)-B1(5,1)*dFtdV*wecs.rb, A1(5,6)+B1(5,1)*dFtdomega ,A1(5,7),0, B1(5,1)*dFtdbeta;...
        A1(6,1:3),A1(6,4)-B1(6,2)*dTrdV,A1(6,5)-B1(6,2)*dTrdV*wecs.rb, A1(6,6)+B1(6,2)*dTrdomega ,A1(6,7),0, B1(6,2)*dTrdbeta;...
        A1(7,1:7), B1(7,3), 0;...
        zeros(1,7),  -1/kappa, 0;...
        zeros(1,8),  -1/tau];
    A7DoF = A9DoF(1:7,1:7);
    
    
    B9DoF = [zeros(3);...
        zeros(1,2),B1(4,1)*dFtdV;...
        zeros(1,2),B1(5,1)*dFtdV;...
        zeros(1,2),B1(6,2)*dTrdV;...
        zeros(1,3);...
        1/kappa, 0,0;...
        0,1/tau,0];
    
    B7DoF = B9DoF(1:7,:);
    B7DoF(:,1) = A9DoF(1:7,8); % Torque Tg as input
    B7DoF(:,2) = A9DoF(1:7,9); % Pitch beta as input
    
    C9DoF = [zeros(1,5),     1     zeros(1,3) % omega_r
        A9DoF(4,:)];% yt_dotdot
    C7DoF = C9DoF(:,1:7);
    
    D9DoF = [zeros(size(C9DoF,1)-1,3); B9DoF(4,:)];
    D7DoF = [zeros(size(C7DoF,1)-1,3); B7DoF(4,:)];
    
    %% Create a state-space model
    if useActuatorStates == 1
        sysLagrange = ss(A9DoF,B9DoF*diag([1000,1,1]),C9DoF,D9DoF); % to give input Tg in kN/m
        Mdl2Str = ' Rot,Twr,Bld+Act ';
    else
        sysLagrange = ss(A7DoF,B7DoF*diag([1000,1,1]),C7DoF,D7DoF); % to give input Tg in kN/m
        Mdl2Str = ' Rot,Twr + Bld';
    end
    
    sysLagrange.Outputname = sysOutputname;
    sysLagrange.Inputname = sysInputname;
    
    sysCell{index} = sysLagrange;
    Torque_g_bar = 0.5*rho*pi*wecs.rb^3*Cqb*Vbar^2/wecs.Ng; % Nm
    u_bar_cell{index} = [Torque_g_bar, beta_bar * pi/180,Vbar]';
    x_bar_cell{index}  = - sysLagrange.A\sysLagrange.B *u_bar_cell{index}; % X(s) = (sI - A)^1 *B *U(s) (xdot = 0 = Ax +Bu)
    y_bar_cell{index}  = sysLagrange.C * x_bar_cell{index} + sysLagrange.D * u_bar_cell{index};
    
    
    %% State space model
    % States omega xdotfa xdotsw xfa xss
    
    % constant for CT: 0.5*rho*A
    kCT = 1/2*rho*pi*wecs.rb^2;
    
    % dotomega = (1/Ng * Tr - Tg)/J
    A11 = 1/wecs.Jr*kCT* wecs.rb*Vbar^2*dCqdlambdab*dlambdadomega;
    A12 = 1/wecs.Jr*kCT*wecs.rb*(Vbar^2*dCqdlambdab*dlambdadxdotfa-2*Vbar*Cqb);
    
    % xdotdotfa = 1/Mt( -Bt * xdotfa - Kt* xfa +Ft)
    A21 = 1/wecs.mt*kCT*Vbar^2*dCtdlambdab*dlambdadomega;
    A22 = 1/wecs.mt*kCT*Vbar^2*dCtdlambdab*dlambdadxdotfa-2*1/wecs.mt*kCT*Vbar*Ctb;
    
    % xdotdotfa = 1/Mt( -Bt * xdotfa - Kt* xfa + 3/(2*H)Tg)
    A31 = 0;%A31 = 1/wecs.mt*kCT*wecs.rb*Vbar^2*dCqdlambdab*dlambdadomega;
    A32 = 0; %A32 = kss/wecs.mt*kCT*wecs.rb*Vbar^2*dCqdlambdab*dlambdadxdotfa-2*kss*1/wecs.mt*kCT*wecs.rb*Vbar*Cqb; % wrong?

    A = [A11 A12 0 0 0;
        A21 (A22-wecs.Bt/wecs.mt) 0 -wecs.Kt/wecs.mt 0;
        A31 A32  -wecs.Bt/wecs.mt 0 -Kss/wecs.mt ;
        [0 1 0 0 0; 0 0 1 0 0]];
    
    % Inputs: Torque [wecs.N/m], pitch angle beta [rad], wind speed V [m/s]
    B11 = -wecs.Ng/wecs.Jr;
    B12 = 1/wecs.Jr*kCT*wecs.rb*Vbar^2*dCqdbetab;
    B13 = 1/wecs.Jr*kCT*wecs.rb*Vbar^2*dCqdlambdab*dlambdadV+2/wecs.Jr*kCT*wecs.rb*Vbar*Cqb;
    
    B22 = 1/wecs.mt*kCT*Vbar^2*dCtdbetab;
    B23 = 1/wecs.mt*kCT*Vbar^2*dCtdlambdab*dlambdadV+2/wecs.mt*kCT*Vbar*Ctb;
        
    B31 = wecs.Ng/wecs.mt *3/(2*wecs.H);
    B32 = 0; %B32 = kss/wecs.mt*kCT*wecs.rb*Vbar^2*dCqdbetab;
    B33 = 0; %B33 = kss/wecs.mt*kCT*wecs.rb*Vbar^2*dCqdlambdab*dlambdadV+kss*2/wecs.mt*kCT*wecs.rb*Vbar*Cqb;

    B = [B11 B12 B13;
        0   B22 B23;
        B31 B32 B33;
        zeros(2,3);];
    
    % Outputs: omega [rad/s], fore-aft dx_dot_f, fore-aft dx_dot_fa [m/s^2]
    C(2,:) = A(2,:);
    % C(3,:) = A(3,:);
    C(1,1) = 1;
    
    D(2,:)=B(2,:);
    %D(3,:)=B(3,:);
    
    sys5DoF = ss(A,B*diag([1000,1,1]),C,D);
    sys5DoF.OutputName = sysOutputname;
    sys5DoF.InputName = sysInputname;
    
    %% Make a Bode plot
    if createBodePlots == 1
        ff = figure(index*100+3 + figNoAdd);%ff.Units='centimeters';
        ff.Position = [520   267   560   531];
        ff.Color = 'white'; ff.Visible = plotVisible;
        
        bodemag(sysLagrange); hold on; grid on;
        % bodemag(ss(A9DoF,B9DoF*diag([1000,1,1]),C9DoF,D9DoF),'k--'); legend('Rot,Twr + Bld', 'Rot,Twr,Bld+Act')
        
        bodemag(sys5DoF,'k--'); hold on; grid on;
        
        for indexModels=1:lenModelNames
            sysm =  sysmCell{indexModels};
            % model = [1 0 0; 0 1 0; 0 0 2*pi/60]*sysm{index}([12 6 38],[8 9 1])*[1000 0 0; 0 1 0; 0 0 1]; %pi/180
            model = [pi/30 0; 0 1]*sysm{index}([38 12],[8 9 1])*[1000 0 0; 0 1 0; 0 0 1]; %pi/180
            model.OutputName = sysOutputname;
            model.InputName = sysInputname;
            bodemag(model,opts); hold on
            grid on
        end
        
        cl = colormap('lines');
        
        legCell = [ {'\color{black}Mdl1: Rot+Twr '};
            {['\color[rgb]{',num2str(cl(1,:)),'}Mdl2:',Mdl2Str]};
            modelsL];
        
        strExtr = '';
        for ilegC = 1: length(modelsL)
            strExtr = [strExtr,...
                '\color[rgb]{',num2str(cl(ilegC+1,:)),'}',' ',modelsL{ilegC}]; %#ok<AGROW> This fast
        end
        
        % legend(legCell,'location','southeast')
        
        titleStr = sprintf('Bode mag. plot for V = %d, pitch = %2.2f, and TSR = %2.4f',Vbar,beta_bar,Lambda_bar);
        titleCell = {titleStr, [legCell{1},legCell{2},strExtr]};
        title(titleCell)
        
        figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),'.png'];
        figFolderStr = fullfile(figFolder,figStr);
        figFolderStrEps = fullfile(figFolder,sprintf('BodeV%d',Vbar));
        % export_fig(figFolderStr,'-dpng','-m2','-nocrop') %ToDo AD Check why
        print(figFolderStr, '-dpng');
        print(figFolderStrEps, '-depsc');
    end
    
end


%% Modal analysis
return; %
% This analysis is left out fro noe
ff = figure(50 + figNoAdd);ff.Units='centimeters'; ff.Position=[2 2 45 27]; ff.Color='white'; %#ok<UNRCH> This analysis is left out fro noe
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

% Antje toDo: Investigate this later
% Ca = [dFtdV, dFtdV*rb, dFtdomega, 0;
%     dFtdV, dFtdV*rb, dFtdomega,0;
%     dTrdV, dTrdV*rb, dTrdomega,0;
%     zeros(1,4)];
% Q1Tilda = [dFtdbeta,dFtdbeta, dTrdbeta, 0]';
%
% A9DoFTest =  [A1(1:3,:), zeros(3,2);
%     [-inv(M)* K, -inv(M) *(Ce + Ca), zeros(4,1), inv(M) *Q1Tilda];
%     zeros(1,7),  -1/kappa, 0;...
%     zeros(1,8),  -1/tau];




