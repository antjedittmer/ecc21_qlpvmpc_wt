%%  Initialize: Add path to data and load LUT for 
clear; close all; clc
addpath(fullfile(fileparts(pwd),'dataIn'));
initModel5MWNREL;

%% Physical and model constants
rho = 1.2250;
rb = 63; % m blade radius
Ng = 97; %Gearbox ratio
Jr = (115926 + 3 * 11.776e6) *1; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
Jgr = Jg * Ng^2; % Inertia of the generator in LSS
Ks = 867637000; %Nm/rad Stiffness of the transmission
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
tau_g = 0.01; % for pitch actuator and airspeed

%% LUT
LUT_V = 4:25; % LUT table rows
LUT_omega = [0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)];
LUT_beta = [zeros(1,8),0.0642    0.1126    0.1493    0.1808    0.2093    0.2356    0.2602...
    0.2833    0.3054    0.3269    0.3477    0.3676    0.3868    0.4056];

%% Lambda, omega and Cq at two wind speeds
% Set to different wind speeds
Vbar = 4; %wind(1);
Vbar2 = Vbar + 1; %wind(2);

% Get rotor speed omegabar at two wind speed, beta = 0;
beta_bar = interp1(LUT_V,LUT_beta,Vbar,'pchip')/pi*180;

% omegaR, lambda at two windspeeds, gradient dlambdadomega at first wind speed 
omegabar = interp1(LUT_V,LUT_omega,Vbar,'pchip');
Lambda_bar = omegabar*rb/Vbar;
dlambdadomega = rb/Vbar;
dlambdadV = - omegabar*rb/Vbar^2;

omegabar2 = interp1(LUT_V,LUT_omega,Vbar2,'pchip');
Lambda_bar2 = omegabar2*rb/Vbar2;

% Cut LUT Rotor_cQ into right format
Rotor_cQ = Rotor_cQ(:,1:length(Rotor_Lambda));
CqPitch0 = Rotor_cQ(Rotor_Pitch == 0,:);    %interp2(Rotor_Lambda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'spline');

%% Compute gradients %%[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
dLambda = mean(diff(Rotor_Lambda));
dPitch = mean(diff(Rotor_Pitch));
[Cqdlambda,CqdbetaDeg] = gradient(Rotor_cQ); %dcQ/dbetaDeg * dbetaDeg/dbetaRad *betaRAd
Cqdbeta = CqdbetaDeg * 180/pi;
Cqdlambda_pitch0 = Cqdlambda(Rotor_Pitch == 0,:);

%% Grid points
dCqdlambdab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdlambda,Lambda_bar,beta_bar,'spline');
dCqdbetab = interp2(Rotor_Lambda,Rotor_Pitch,Cqdbeta,Lambda_bar,beta_bar,'spline');
Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'spline');
Cqb2 = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar2,beta_bar,'spline');

%% Caltulate Tr=0.5*pi*rb^3*Cq*v^2
dTrdlambda = 0.5*rho*pi*rb^3* dCqdlambdab *Vbar^2;
dTrdomega = 0.5*rho*pi*rb^3*Vbar^2*(dCqdlambdab*dlambdadomega);% Bianchi: br = dTr/domega
dTrdbeta = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;% Krb = dTr/dbeta
dTrdV =  2*0.5*rho*pi*rb^3*Vbar*Cqb + ...
    0.5*rho*pi*rb^3*Vbar^2*dCqdlambdab*dlambdadV;

Tr =  0.5*rho*pi*rb^3 * Vbar^2 * Cqb;
Tr2 =  0.5*rho*pi*rb^3 * Vbar2^2 * Cqb2;

Trdiff0 = dTrdlambda * (Lambda_bar2 - Lambda_bar);
TrDiffOmega  = dTrdomega * (omegabar2 - omegabar) ;
TrDiffV = dTrdV *(Vbar2 - Vbar);

dTr = TrDiffOmega + TrDiffV;

%% Display results torque
% Tr(lambda2) = Tr(lambda1) + [dTr(lambda)/dlambda| lambda1] * (lambda
fprintf('Difference torque: Tr(wind2) - Tr(wind): %2.3e,  Sum gradient terms: %2.3e \n', Tr2 - Tr, dTr)
fprintf('Gradient term dTrdomega * (omegabar2 - omegabar): %2.3e \n', dTrdomega * (omegabar2 - omegabar))
fprintf('Gradient term dTrdV *(Vbar2 - Vbar): %2.3e \n\n', dTrdV *(Vbar2 - Vbar))

%% Display results torque
fprintf('Difference torque coefficient: Cq(wind2) - Cq(wind): %2.3e, with Cq(wind): %2.3e, Cq(wind2): %2.3e \n', Cqb2 - Cqb, Cqb, Cqb2)
fprintf('Ratio (Cq(wind2) - Cq(wind))/(lambda2 -lambda1) %2.3e \n', (Cqb2 - Cqb)/(Lambda_bar2 -  Lambda_bar))
fprintf('Gradient at dCqdlambdab %2.3e \n', dCqdlambdab)
fprintf('Cq(wind2), Cq(wind1) + dCqdlambdab * (lambda2 - lambda1): %2.3e \n', Cqb2, Cqb + dCqdlambdab* (Lambda_bar2 -  Lambda_bar))
fprintf('Cq(omega2), Cq(omega1) + dCqdlambdab * (lambda2 - lambda1): %2.3e \n', Cqb2, Cqb + dCqdlambdab*dlambdadomega*(omegabar2 - omegabar))

%% Plot Torque coefficient vs. lambda
figure(1); 
ax(1) = subplot(2,1,1);
plot(Rotor_Lambda,CqPitch0, Lambda_bar,Cqb,'k*',Lambda_bar2,Cqb2,'ro')
axis tight;  grid on;
xlabel('TSR \lambda [-]'); ylabel('Cq_{\beta 0}(\lambda)')
title('Torque coefficient Cq_{\beta 0}(\lambda) (Pitch angle \beta = 0)')
legend('Cq at Pitch \beta = 0', 'Cq at wind V_1', 'Cq at wind V_2');
ax(2) = subplot(2,1,2);
plot(Rotor_Lambda,Cqdlambda_pitch0, Rotor_Lambda(2:end), CqPitch0(2:end) - CqPitch0(1:end-1),'--',...
    Lambda_bar,dCqdlambdab,'k*'); axis tight;  grid on;
legend('\delta Cq(\lambda)/ \delta \lambda at Pitch \beta = 0', ...
    'Cq(2:end) - Cq(1:end-1)', '\delta Cq(\lambda)/ \delta \lambda at wind V_1');
xlabel('TSR \lambda[-]'); ylabel('\delta Cq_{\beta 0}(\lambda)/ \delta \lambda')
linkaxes(ax,'x');

