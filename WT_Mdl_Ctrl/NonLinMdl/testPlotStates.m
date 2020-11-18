close all;
%load('OutDataTestOL_LPV_windChange_20sec', 'xDot06','xkLPV','xk06', 'u','uDot','A','B','wind');
%load('OutDataTestOL_LPV_windChange_200sec', 'xDot06','xkLPV','xk06', 'u','uDot','A','B','wind');
load('OutDataTestOL_LPV_windChange_400sec', 'xDot06','xkLPV','xk06', 'u','uDot','A','B','wind');


xCell = {'\omega_{g} [rad/s]','T_g [Nm]','\Delta \Theta_s [rad]','\Delta \omega_r [rad/s]',...
    '\Delta \omega_{gr} [rad/s ]','\Delta T_g [Nm]'};

figure(1);

ax1 = nan(8,1);
%DT = 2* DT;
time = 0: DT : (size(xk06,1)-1)*DT;
ax1(1) = subplot(8,1,1);
plot(time, wind(1:length(time)));
ylabel('wind [m/s]');
axis tight; grid on;
title('Comparison NonLin/ LPV model')

ax1(2) = subplot(8,1,2);
plot(time, u(:,2)*180/pi);
ylabel('pitch \beta [deg]');
axis tight; grid on;


for idx = 1:6
    ax1(idx+2) = subplot(8,1,idx+2);
    plot(time,xk06(:,idx),time, xkLPV(:,idx),'--' );
    ylabel(xCell{idx})
    axis tight; grid on;
end
xlabel('time [s]')
legend('NonLin','LPV');
linkaxes(ax1,'x')

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.5]);

% rNew = [r;zeros(4,1)];
%figure; for idx = 1:7, subplot(7,1,idx); myXk = Xk(idx:7:end); plot(myXk); hold on; plot(ref(idx:7:end),'r--');ylabel(num2str(myXk(end))); end
figure; for idx = 1:8, subplot(8,1,idx); myXk = Xk(idx:8:end); plot(myXk); hold on; plot(ref(idx:8:end),'r--');ylabel(num2str(myXk(end))); end

figure; for idx = 1:2, subplot(2,1,idx); myXk = UDot(idx:2:end); plot(myXk); ylabel({num2str(myXk(1)),num2str(myXk(end))}); end

%
%load('OutDataTestOL_LPV', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B');
%load('OutDataTestOL_LPV_windConst5', 'xDot06','OutDataTestOL_LPV','xk06', 'u','uDot','A','B');

% rho = 1.2250;
% rb = 63; % m blade radius
% Lambda_bar = omegabar*rb/Vbar;
% Ng = 97; %Gearbox ratio
% Jr = 115926 + 3 * 11.776e6; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades)
% Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
% Jgr = Jg * Ng^2; % Inertia of the generator in LSS
% Ks = 867637000; %Nm/rad Stiffness of the transmission
% Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
% tau_g = 0.01; % for pitch actuator and airspeed
% 
% %% xk06: States nonlinear model [omegaG,Tg,ThetasDot,omegaRdot,omegaGdot,TgDot]
% omegaGRdotVec_kplus1 = xk06(2:end,5);
% ThetaOmegaROmegaGrTg_k =  xk06(1:end-1,3:6);
% vecNonLin = ThetaOmegaROmegaGrTg_k \ (omegaGRdotVec_kplus1* Jgr*Ng);test
% 
% ThetadotVec_kplus1 = xk06(2:end,3);
% vecNonLinTheta = ThetaOmegaROmegaGrTg_k \ (ThetadotVec_kplus1);
% 
% mat_kplus1 = xk06(2:end,3:6);
% matrixNonLin = ThetaOmegaROmegaGrTg_k \ mat_kplus1;
% 
% 
% omegaGRdotVec_kplus1Lin = OutDataTestOL_LPV(2:end,5);
% ThetaOmegaROmegaGrTg_kLin =  OutDataTestOL_LPV(1:end-1,3:6);
% vecLin = ThetaOmegaROmegaGrTg_kLin \ (omegaGRdotVec_kplus1Lin* Jgr *Ng);
% 
% 
% % vec sys matrix: Ks/Jgr, Bs/Jgr, -Bs/Jgr, Ng/Jgr;
% fprintf('Ks: Phys. value: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  Ks, vecNonLin(1), vecLin(1))
% fprintf('Bs: Phys. value: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  Bs, vecNonLin(2), vecLin(2))
% fprintf('-Bs: Phys. value: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  -Bs, vecNonLin(3), vecLin(3)/Ng)
% fprintf('Ng: Phys. value: %2.3e, NonLin  %2.3e, Lin %2.3e \n', Ng, vecNonLin(4), vecLin(4))