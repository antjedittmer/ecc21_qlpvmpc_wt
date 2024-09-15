clear; clc; close all
% addpath('C:\Program Files\export_fig'); %C:\Users\bs\Desktop\export_fig-master\export_fig-master;
% addpath C:\Werk_Data\Software\hline_vline;
% load 'NREL5MW_linearised_reduced'
addpath(fullfile(fileparts(pwd),'dataIn'))
load('NREL5MW_linearised')
load 'NREL5MW_CPdata_new3'
Lin_points = Lin;

%% Plot parameters
IPCplots = 0;
plotVisible = 'on';
sysOutputname = [{'omega [rad/s]'}; {'xdotdot_{fa} [m/s^2]'}; {'xdotdot_{ss} [m/s^2]'}];
sysInputname = [{'Torque [kNm]'}; {'\beta [rad]'}; {'V [m/s]'}];

%% Create output folder
figFolder = 'FigureNewLayoutTUDelft_201014'; %'Figures';
if ~isfolder(figFolder)
    mkdir(figFolder);
end

%% Define parameters
J = 4.2e7;% 3.5444E+07
R = 63;
rho = 1.225;
G = 97;
Effic = 0.944;
vecLambda = 1:150;
Rotor_Lambda = Rotor_Lamda(vecLambda);
Rotor_Pitch = Rotor_Pitch * pi/180; % Rotor Pitch deg ->rad
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'cubic');
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'cubic');
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'cubic');

%% Compute gradients
dX = mean(diff(Rotor_Lambda));
dY = mean(diff(Rotor_Pitch));
[Cqdlambda,Cqdbeta] = gradient(Cq,dX,dY);
[Cpdlambda,Cpdbeta] = gradient(Cp,dX,dY);
[Ctdlambda,Ctdbeta] = gradient(Ct,dX,dY);
Clambda=Cq./Xq.^2;
[Clambdadlambda,Clambdadpitch] = gradient(Clambda,dX,dY);

%% Initialize Wnn ZZ PP matrices
sizeA = 5; % number of states, ToDo AD: Is there an automatic way to determine this?
Wnn = NaN(sizeA,length(Lin.V));
ZZ = NaN(sizeA,length(Lin.V));
PP = NaN(sizeA,length(Lin.V));

%% Parameters for tower model
kfa = 1; %should this be a plus or a min?
kss = 0.1; % correction factor -> this can be removed with 2/(3H) ToDo
mt = 4.38e5; % mt+N*mb: 400680
ct = 2*6421;
Kt = 1.846e6;

%% Options for Bode plots
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz'; opts.MagUnits = 'Abs'; opts.MagScale ='log'; opts.XLim ={[1e-2 1e1]};

%% Initialize name of models and load
idxModelNames = 1:2; %1:5

modelNames =[{'NREL5MW_linearised'}; {'NREL5MW_linearised_Gen_Fa_Flap1'};
    {'NREL5MW_linearised_Gen_Fa_SS_Flap1'};
    {'NREL5MW_linearised_Gen_DriveTrain_Fa_SS_Flap1'}; {'NREL5MW_linearised_reduced'} ];

modelsL=[{'FAST (30 states)'};
    {'FAST Gen+FA+1^{st} flap mode (10 states)'};
    {'FAST Gen+FA+SS+1^{st} flap mode (12 states)'};
    {'FAST Gen+drive train+FA+SS+1^{st} flap mode (14 states)'};
    {'FAST Gen+FA (4 states)'}];

modelNames = modelNames(idxModelNames);
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
sysCell = cell(lenLinV,1);

%% Loop over all air speeds in LinV
for index =  length(Lin.V) %[1, length(Lin.V)]%
    
    %% Variables dependent on wind speed V
    Vbar = Lin_points.V(index);
    omegabar = Lin_points.RSpeed(index);
    beta_bar = Lin_points.Pitch(index); %/pi*180;
    
    Lambda_bar= omegabar*R/Vbar;
    dlambdadomega = R/Vbar;                 %s
    dlambdadV = -omegabar*R/Vbar^2;         %s/m
    dlambdadxdotfa = -kfa*omegabar*R/Vbar^2; %s/m   xdotfa m/s
    
    %% Grid points
    dCqdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
    dCtdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdlambda,Lambda_bar,beta_bar,'cubic');
    dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
    dCtdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic');
    Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
    Ctb = interp2(Rotor_Lambda,Rotor_Pitch,Ct,Lambda_bar,beta_bar,'cubic');
    
    %% State space model
    % dotomega=(Tr-GTg)/J
    % States omega xdotfa xdotss xfa xss
    kCT = 1/2*rho*pi*R^2;
    A11 = 1/J*kCT*R*Vbar^2*dCqdlambdab*dlambdadomega;
    A12 = 1/J*kCT*R*Vbar^2*dCqdlambdab*dlambdadxdotfa-2*kfa/J*kCT*R*Vbar*Cqb;
    
    A21 = 1/mt*kCT*Vbar^2*dCtdlambdab*dlambdadomega;
    A22 = 1/mt*kCT*Vbar^2*dCtdlambdab*dlambdadxdotfa-2*kfa/mt*kCT*Vbar*Ctb;
    
    A31 = kss/mt*kCT*R*Vbar^2*dCqdlambdab*dlambdadomega;
    A32 = kss/mt*kCT*R*Vbar^2*dCqdlambdab*dlambdadxdotfa-2*kss*kfa/mt*kCT*R*Vbar*Cqb; % wrong
    
    A = [A11 A12 0 0 0;
        A21 (A22-ct/mt) 0 -Kt/mt 0;
        A31 A32  -ct/mt 0 -Kt/mt ;
        [0 1 0 0 0; 0 0 1 0 0]];
    
    % Inputs: Torque [N/m], pitch angle beta [rad], wind speed V [m/s]
    B11 = -G/J;
    B12 = 1/J*kCT*R*Vbar^2*dCqdbetab;
    B13 = 1/J*kCT*R*Vbar^2*dCqdlambdab*dlambdadV+2/J*kCT*R*Vbar*Cqb;
    B22 = 1/mt*kCT*Vbar^2*dCtdbetab;
    B23 = 1/mt*kCT*Vbar^2*dCtdlambdab*dlambdadV+2/mt*kCT*Vbar*Ctb;
    
    B32 = kss/mt*kCT*R*Vbar^2*dCqdbetab;
    B33 = kss/mt*kCT*R*Vbar^2*dCqdlambdab*dlambdadV+kss*2/mt*kCT*R*Vbar*Cqb;
    
    B=[B11 B12 B13;
        0   B22 B23;
        0   B32 B33;
        zeros(2,3);];
    
    % Outputs: omega [rad/s], fore-aft dx_dot_f, fore-aft dx_dot_fa [m/s^2]
    C(2,:) = A(2,:);
    C(3,:) = A(3,:);
    C(1,1) = 1;
    
    D(2,:)=B(2,:);
    D(3,:)=B(3,:);
    
    %% Make a Bode plot
    sys = ss(A,B*diag([1000,1,1]),C,D); 
    
    ff = figure(index*100+2); %ff.Units='centimeters'; ff.Position=[2 2 45 27];
    ff.Position = [520   267   600   600];
    ff.Color='white'; ff.Visible = plotVisible;
    sys.Outputname = sysOutputname;
    sys.Inputname = sysInputname;
    sysCell{index} = sys;
    
    bodemag(sys,'k'); hold on
    
    for indexModels=1:lenModelNames
        sysm =  sysmCell{indexModels};
        model = [2*pi/60 0 0; 0 1 0; 0 0 1]*sysm{index}([38 12 13],[8 9 1])*...
             [1000 0 0; 0 1 0; 0 0 1]; %[1000 0 0; 0 (1/180*pi) 0; 0 0 1];
        model.OutputName = sysOutputname;
        model.InputName = sysInputname;
        bodemag(model,opts); hold on
        grid on
    end
    
    [strExtr,legCell] = createLegendCell(modelsL);
    legCell{1} = 'TU Delft 5 states';
    %legend(legCell,'location','southeast')
    titleStr = sprintf('Bode mag. plot for V = %d, pitch = %2.2f, and TSR = %2.4f',Vbar,beta_bar,Lambda_bar);
    titleCell = {titleStr,[legCell{1},strExtr]};
    title(titleCell)
    figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),'.png'];
    figFolderStr = fullfile(figFolder,figStr);
    % export_fig(figFolderStr,'-dpng','-m2','-nocrop') %ToDo AD Check why
    print(figFolderStr, '-dpng');
    
    %% Make a Bode plot IPC (no MBC)
    if IPCplots == 1
        ff = figure(index*100+2); %ff.Units='centimeters'; ff.Position=[2 2 45 27];
        ff.Position = [520   267   560   531];
        ff.Color='white';ff.Visible=plotVisible;
        for indexModels = 1:lenModelNames
            sysm = sysmCell{indexModels};
            bodemag(sysm{index}([21 22 23],[4 5 6]),opts); hold on
        end
        grid on
        
        %legend(legCell,'location','southeast')
        titleStr = strrep(titleStr,'Bode mag. plot for V','Bode mag. for IPC (no MBC) plot for V');
        titleCell = {titleStr,[legCell{1},strExtr]};
        title(titleCell)
        
        figStr = strrep(figStr,'Figure_Bode V =','Figure_Bode_for_IPC_nMBC_V =');
        figFolderStr = fullfile(figFolder,figStr);
        %export_fig(figFolderStr,'-dpng','-m2','-nocrop')
        print(figFolderStr, '-dpng');
        
        %% Make a Bode plot IPC MBC
        ff = figure(index*100+3);%ff.Units='centimeters'; ff.Position=[2 2 45 27]; 
        ff.Position = [520   267   560   531];
        ff.Color='white';
        ff.Visible = plotVisible;
        for indexModels = 2:lenModelNames
            sysm_MBC = sysmCell_MBC{indexModels};
            bodemag(sysm_MBC{index}([21 22 23],[4 5 6]),opts); hold on
        end
        grid on
        
        %legend(legCell,'location','southeast')
        titleStr = strrep(titleStr,'Bode mag. plot for V','Bode mag. for IPC (no MBC) plot for V');
        titleCell = {titleStr,[legCell{1},strExtr]};
        title(titleCell)
        
        legend(modelsL{2:end},'location','southeast')
        titleStr = strrep(titleStr,'Bode mag. for IPC (no MBC) plot for V','Bode mag. for IPC MBC plot for V');
        title(titleStr);
        
        figStr = strrep(figStr,'Figure_Bode_for_IPC_nMBC_V =','Figure_Bode_for_IPC_MBC_V = ');
        figFolderStr = fullfile(figFolder,figStr);
        %export_fig(figFolderStr,'-dpng','-m2','-nocrop')
        print(figFolderStr, '-dpng');
    end
    
    %% Look at pole locations
    ff = figure(index*100+4);
    %ff.Units='centimeters'; ff.Position=[2 2 45 27]; 
    ff.Color='white'; ff.Visible = plotVisible;
    poleMarkerSize = 10;
    [Wnn(:,index),ZZ(:,index),PP(:,index)]=damp(sys);
    plot(real(PP),imag(PP),'k+','MarkerSize',poleMarkerSize); hold on
    symbolStr = {'o','s'};
    cl = colormap('lines');
    for indexModels = 1:lenModelNames
        sysm = sysmCell{indexModels};
        sysm_MBC = sysmCell_MBC{indexModels};
        [Wn{index,indexModels},Z{index,indexModels},P{index,indexModels}]=damp(sysm{index});
        [Wnmbc{index,indexModels},Zmbc{index,indexModels},Pmbc{index,indexModels}]=damp(sysm_MBC{index});
        plot(real(P{index,indexModels}),imag(P{index,indexModels}),...
            symbolStr{indexModels},...
            'MarkerSize',poleMarkerSize,'Color',cl(indexModels,:),...
            'LineStyle','none'); hold on
    end
    grid on;
    % legend([{'Simple model Drive-Train, FA, SS (5 states)'};modelsL],'location','southeast')
    titleStr = strrep(titleStr,'Bode mag. plot for V =','Pole locations. for IPC  plot for V');
    titleStr = strrep(titleStr,'Bode mag. for IPC (no MBC) plot for V','Pole locations. for IPC  plot for V');
    titleCell = {titleStr,[legCell{1},strExtr]};  
    title(titleCell);
    
    
    figStr = strrep(figStr,'Figure_Bode V =', 'Figure_Poles_V =');
    figStr = strrep(figStr,'Figure_Bode_for_IPC_nMBC_V =','Figure_Poles_V = ');
    figFolderStr = fullfile(figFolder,figStr);
    %export_fig(figFolderStr,'-dpng','-m2','-nocrop')
    print(figFolderStr, '-dpng');
    
end

% save sysCell to outputfolder
save(fullfile(figFolder,'sysCellLin_5states.mat'),'sysCell');

% if this is empty return to prevent crash
if isempty(Zmbc{1,4})
    return;
end


ff = figure(50); % ff.Units='centimeters'; ff.Position=[2 2 45 27];
ff.Color='white';
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

% gcf1 = figure(1);
% bodemag(sys,'k'); hold on;
%
% tfSYS = tf(sys);
% tf_torque2omega = tfSYS(1,1);
% tf_torqueBetaV2omega = tfSYS(1,:);
% tf_torqueBetaV2xdotdotFA = tfSYS(2,:);
% tf_V2omegaBetadV = tfSYS(3,:);
%
% figure(2);
%
% % Omega: 1st order dynamic %zpk_torque2omega = zpk(tf_torque2omega);
% for idxPz = 1:3
%     pzmap(tf_torqueBetaV2omega(1,idxPz));
%     hold on,
% end
% legend(sysInputname)
% poleOmega = pole(tf_torqueBetaV2omega);
% minPole = min(poleOmega);
% dcgain(tf_torqueBetaV2omega);
%
% figure(4);
% for idxPz = 1:3
%     pzmap(tf_torqueBetaV2xdotdotFA(1,idxPz));
%     hold on,
% end
%
% ATest = [1,0; -Cfa/Mfa -Kfa/Mfa];  BTest = [0,1/Mfa]';
% CTest = ATest(2,:); DTest = BTest(2);
% ssTest = ss(ATest,BTest,CTest, DTest);
