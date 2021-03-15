function Linear_Models_Lagrange_Comparison(speedVec,figFolder, useActuatorStates,figNoAdd,createBodePlots,plotVisible)

%% Set defaults
if nargin == 0
    speedVec = [1,22]; %[1,22]; % variations of wind speed
    figFolder = 'FigureLagrange_210308_comp';
    useActuatorStates = 0;
    figNoAdd = 0;
    close all;
    createBodePlots = 1;
    plotVisible = 'on';
end

%% Initialize: Set path
addpath(fullfile(fileparts(pwd),'dataIn'));
load('NREL5MW_CPdata','Rotor_Lamda', 'Rotor_Pitch', 'Rotor_cQ', 'Rotor_cT');

%% Initialize name of models and load into workspace
idxModelNames = 1; %1:2; %1:5
modelNames ={'NREL5MW_linearised_4to25'};
modelNames = modelNames(idxModelNames);

modelsL = {'FAST 30 states'};
modelsL = modelsL(idxModelNames);

lenModelNames = length(modelNames);
sysmCell = cell(lenModelNames,1);

idxModels = 1;

load(modelNames{idxModels},'sysm','Lin'); %'sysm_MBC');
sysmCell{idxModels} = sysm;
for idxModels = 2 : lenModelNames
    load(modelNames{idxModels},'sysm');
    sysmCell{idxModels} = sysm;
end


%% Plot parameters
% yt is the axial displacement of the nacelle, $\zeta$ is the angular blade displacement
% out of the plane of rotation and $\omega_r$ and $\omega_r$ are the angular positions of the
% X =[yt zeta thetha_r thetha_gr yt_dot zeta_dot omega_r omega_gr ]'
% uact=[Tg_ref beta_ref V ]'
% y = [ yt_dot zeta_dot omega_r omega_gr]'
% u = [beta_ref wg_ref]' d=v
% y = [yt_dot zeta_dot omega_r omega_z]'

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
[wecs, M, Ce, K, Q, L, rho, tau, kappa, Rotor_Lambda, Rotor_Pitch, Cq, Ct, Cqdlambda,Cqdbeta, Ctdlambda,Ctdbeta] = initModel5MWNREL(0, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figFolder);

%% Constant Matrices 
A1 = [zeros(4), L;...
    -M\K, -M\Ce];

B1 = [zeros(4); M\Q];

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
for index =  speedVec %[1,9,10,length(Lin.V)]%1 : length(Lin.V) %[1,length(Lin.V)] % [1,2,length(Lin.V)-1,length(Lin.V)]
    
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
    dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic'); %* 180/pi;
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
    % Models from Bianchi do not include sidevards motions
    % X =[yt zeta Thetha_s yt_dot zeta_dot omega_r omega_gr beta torque]'
    % y=[ yt_ddot rb*zeta_dot omega_r]' -> y=[ yt_ddot rb*zeta_dot omega_r]'
    % uact=[Tg_ref beta_ref V]'
    % Ve = V - yt_dot [m/s] - rb [m]  * zeta_dot [1/s]
    % n
    
    ksw = 3/(2* wecs.H);
    
    A9DoF = [A1(1:4,:) ,zeros(4,2); % states 1 - 4: x_FA,zeta,x_sw, Theta_s
        [A1(5,1:4), ... % state 5: xdot_FA (from  states 1 - 4: x_FA,zeta,x_sw, Theta_s)
        A1(5,5) - B1(5,1)*dFtdV,... % state 5-> 5: xdot_FA <- xdot_FA
        A1(5,6) - B1(5,1)*dFtdV*wecs.rb, ...% state 6-> 5: xdot_FA <- zeta
        A1(5,7), ... % state 7 > 5: xdot_FA <- xdot_sw
        A1(5,8) + B1(5,1)*dFtdomega, ... % state 8 > 5: xdot_FA <- omega_r
        A1(5,9), ...  % state 9 > 5: xdot_FA <- xdot_sw % state 10, 11 > 5: xdot_FA <- Tg, beta
        0 , B1(5,1)*dFtdbeta];
        [A1(6,1:4), ...% state 6: zeta (from  states 1 - 4: x_FA,zeta,x_sw, Theta_s)
        A1(6,5) - B1(6,1)*dFtdV,... % state 5> 6: zeta<- xdot_FA
        A1(6,6) - B1(6,1)*dFtdV*wecs.rb,... % state 5> 6: zeta<- xdot_FA
        A1(6,7), ... % zeta<- xdot_FA
        A1(6,8) + B1(6,1)*dFtdomega, ... % state 5> 6: zeta<- xdot_FA
        A1(6,9), ... % state 5> 6: zeta<- xdot_FA
        0 ,B1(6,1)*dFtdbeta];
        [A1(7,1:4), ... % state 6: x_sw (from  states 1 - 4: x_FA,zeta,x_sw, Theta_s)
        A1(7,5) - ksw * B1(7,2)*dTrdV,...
        A1(7,6) - ksw * B1(7,2)*dTrdV*wecs.rb,...
        A1(7,7), ...
        A1(7,8) + ksw * B1(7,2)*dTrdomega, ...
        A1(7,9),...
        0 ,ksw * B1(7,2)*dTrdbeta];
        [A1(8,1:4), ...
        A1(8,5) - B1(8,3)*dTrdV,...
        A1(8,6) - B1(8,3)*dTrdV*wecs.rb, ....
        A1(8,7), ...
        A1(8,8) + B1(8,3)*dTrdomega, ...
        A1(8,9), ...
        B1(8,1),B1(8,3)*dTrdbeta];
        [A1(9,1:6), A1(9,7) , A1(9,8:9) ,B1(9,4) ,0];
        zeros(1,9),  -1/kappa, 0;...
        zeros(1,10),  -1/tau];
    
    A7DoF = A9DoF(1:9,1:9);
       
    B9DoF = [zeros(4,3);... % states 1 - 4: x_FA,zeta,x_sw, Theta_s
        zeros(1,2),B1(5,1)*dFtdV;... % state 5: xdot_FA
        zeros(1,2),B1(6,1)*dFtdV;...% state 6: zetadot
        zeros(1,2),B1(7,2)*ksw *dTrdV;...% state 7: xdot_SW
        zeros(1,2),B1(8,3)*dTrdV;... state 8: omega_r
        zeros(1,3);... state 9: omega_gr
        1/kappa, 0,0;...
        0,1/tau,0];
    
    B7DoF = B9DoF(1:9,:);
    B7DoF(:,1) = A9DoF(1:9,8+2); % Torque Tg as input
    B7DoF(:,2) = A9DoF(1:9,9+2); % Pitch beta as input
    
    C9DoF = [zeros(1,5+2),     1     zeros(1,3) % state 8: omega_r
        A9DoF(5,:)]; % state 5: xdot_FA: output xdotdot_FA
    C7DoF = C9DoF(:,1:9);
    
    D9DoF = [zeros(size(C9DoF,1)-1,3); B9DoF(5,:)];
    D7DoF = [zeros(size(C7DoF,1)-1,3); B7DoF(5,:)];
    
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
    
    Cout = [C9DoF;
        [zeros(1,7+2) 1 zeros(1,1)];  % beta
        [zeros(1,8+2) 1]];% torque
    Dout = [D9DoF; zeros(2,3)];
    sysOut = ss(A9DoF,B9DoF*diag([1000,1,1]),Cout,Dout); % to give input Tg in kN/m
    
    sysOut.Outputname = [sysOutputname;{'beta'};{'Tg'}];
    sysOut.Inputname = sysInputname;
    
    
    sysCell{index} = sysOut;
    Torque_g_bar = 0.5*rho*pi*wecs.rb^3*Cqb*Vbar^2/wecs.Ng; % Nm
    u_bar_cell{index} = [Torque_g_bar, beta_bar * pi/180,Vbar]';
    x_bar_cell{index}  = - sysLagrange.A\sysLagrange.B *u_bar_cell{index}; % X(s) = (sI - A)1 *B *U(s) (xdot = 0 = Ax +Bu)
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
        A31 A32  -wecs.Bt/wecs.mt 0 -wecs.Kt/wecs.mt ;
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
    
    D(2,:)= B(2,:);
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
        
%         figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),'.png'];
%         figFolderStr = fullfile(figFolder,figStr);
%         figFolderStrEps = fullfile(figFolder,sprintf('BodeV%d',Vbar));
%         % export_fig(figFolderStr,'-dpng','-m2','-nocrop') %ToDo AD Check why
%         print(figFolderStr, '-dpng');
%         print(figFolderStrEps, '-depsc');
        figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),'.png'];
        figFolderStr = fullfile(figFolder,figStr);
        figFolderStrEps = fullfile(figFolder,sprintf('BodeV%d',Vbar));
        posaxes =  get(gcf,'Position');
        set(gcf,'Position',[posaxes(1:3),posaxes(4)*0.75]);
        title(sprintf('V = %d: %s', Vbar, titleCell{2}))
        % export_fig(figFolderStr,'-dpng','-m2','-nocrop') %ToDo AD Check why
        print(figFolderStr, '-dpng');
        print(figFolderStrEps, '-depsc');

    end
    
end

save('sysCell.mat','sysCell');

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
