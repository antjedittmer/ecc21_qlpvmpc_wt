close all;
addpath(genpath(pwd))
rhoP = 18; %[V,omega, beta];

%% The same for the linear wind turbine model
LinPointsMFile;
useFASTmdl = 0;

% Get time and input vectors to simulate with lsim
t = 0:0.02:100;
t1L = round(length(t)/2);
V = [rhoP*ones(1,t1L),(rhoP+2) * ones(1,length(t) - t1L) ];
d = V + (rand(size(t)) - mean(rand(size(t)))); 

u1 = nan(size(d));
u2 = nan(size(d));
for idx = 1: length(V)
    u1(idx) = Lin_points.Torque(Lin_points.V == V(idx))*10^3; % * ones(size(d))*10^3;
    u2(idx) = Lin_points.Pitch(Lin_points.V == V(idx)); % * ones(size(d));
end

% The input and output vectors are:
% u = [Tg_ref beta_ref V]': sysOut.InputName' {'Torque [kNm]'}    {'\beta [rad]'}    {'V [m/s]'}
% y = [yt_dot zeta_dot omega_r omega_gr]'
% {'omega [rad/s]'}    {'ydotdot_{fa} [m/s^2]'}    {'beta'}    {'Tg'}
% Mdl2 is a 9th order model with state vector:
% x = [yt zeta xt theta_s yt_dot zeta_dot xt_dot omega_r omega_gr]'
% zeta the angular blade out-of-plane displacement and theta_s is the slip
% angle. omega_gr is the generator speed in the low speed shaft reference
% frame. Input u and output y are defined as above

if useFASTmdl == 1
    
    %load('sysOut.mat','sysOut');
    
    sysOut = createLinearModels(rhoP-3,1);
    
    sysOut.StateName = {'yt [m]', 'zeta [m]', 'xt [m]', 'theta_s [rad]', ...
        'yt_dot [m/s]', 'zeta_dot [m/s]', 'xt_dot [m/s]',...
        'omega_r [rad/s]', 'omega_gr [rad/s]','Tg [kNm]','beta [rad]'};
    
    A = sysOut.A; B = sysOut.B(:,1:2); C = sysOut.C(1:4,:); E = sysOut.B(:,3);
      
    % Simulate with lsim
    t = 0:0.02:100;
    t1L = round(length(t)/2);
    
    rhoP1 = Lin_points.RSpeed(Lin_points.V == rhoP); % %omega
    rhoP2 = Lin_points.Pitch(Lin_points.V == rhoP);
   
    u1 = nan(size(d));
    u2 = nan(size(d));
    for idx = 1: length(V)
        u1(idx) = Lin_points.Torque(Lin_points.V == V(idx))*10^3; % * ones(size(d))*10^3;
        u2(idx) = Lin_points.Pitch(Lin_points.V == V(idx)); % * ones(size(d));
    end
    
    %{'theta_s [rad]','omega_r','omega_gr','Tg','beta'};
    x0 = [zeros(1,7),rhoP1,rhoP1,u1(1),u2(1),V(1)]';
    %{'omega_r','omega_gr','Tg','beta','V'};
    y0 = x0(2:end);
    
    
else
    
    LUTMFile;
    MdlMFile;
    clear sysOut;
    
    rhoP1 = Lin_points.RSpeed(Lin_points.V == rhoP); % * ones(size(d))*10^3;
    rhoP2 = Lin_points.Pitch(Lin_points.V == rhoP);
    rhoP123 = [rhoP1,rhoP2, rhoP]
    x0 = [0,rhoP1,rhoP1,u1(1),u2(1),V(1)]';
    y0 = x0(2:5);
    
    [Av,Bv] = ABqL_WT_beta(rhoP123,lambda,beta,Cq);
   
%   testX = -inv(Av) * Bv * u0'; testX(1:3)  
    sysOut.StateName = {'theta_s [rad]','omega_r [rad/s]','omega_{gr} [rad/s]','Tg [Nm]','beta [rad]'};
    sysOut.OutputName = [sysOut.StateName(2:5),{'V [m/s]'}];
    
    A = Av; B = Bv(:,2:3);
    nx = size(Av,2);
    C = [zeros(nx-1,1),eye(nx-1)]; %
    E = Bv(:,1);
    
%     sysTest = ss(Av,Bv,C,zeros(size(C,1),size(Bv,2)));
%     u0 = [V(1),u1(1),u2(1)];
%      Atest = - Bv*u0' / x0(1:5);
   
end



%% Check rank
fprintf('WECS: rank(C*E) is %d; rank(E) is %d\n', rank(C*E),rank(E));
nx = size(A,1);
% nd = size(E,2); %nu = 0;
ny = size(C,1);

% Also create an additional state for the disturbance
A1 =[[A, zeros(nx,1)]; [zeros(1,nx),-1]];
C1 = [C, zeros(ny,1)]; %[zeros(1,nx),1]];
E1 = [E;1]; %xd_dot = -xd + Ed
B1 = [B; zeros(1,2)];
fprintf('rank(C*E) is %d; rank(E) is %d\n', rank(C1*E1),rank(E1));

A = A1;
C = C1;
E = E1;
B = B1;

%% Calculate observer matrices
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
sys = ss(A,[E,B],C,zeros(ny,nd+2));

% Est with no control input: zdot = Fz + Ky, xhat = z + Hy
sysEst = ss(F,[T*B,K],eye(nx),[zeros(size(T*B)), H]);

in1 = [d;u1;u2];
diffIn = diff(in1')';
%[y,~,x] = lsim(sys,[d;u1;u2],t,x0);
[dy,~,dx] = lsim(sys,diffIn,t(1:end-1));
ytemp = y0'+cumsum(dy);
y = [ytemp;ytemp(end,:)];
xtemp = x0' + cumsum(dx);
x = [xtemp;xtemp(end,:)];

%z0 = T*x0; %(eye(size(H*C))- H*C); 
%[xest,~,z] = lsim(sysEst,[u1;u2;y'],t,z0);
inEst = [u1;u2;y'];
diffInEnst = diff(inEst')';
[dXest,~,z] = lsim(sysEst,diffInEnst,t(1:end-1));
xtempEst = x0' + cumsum(dXest);
xest = [xtempEst;xtempEst(end,:)];

figure(1)
ny = size(y,2);
for idx = 1:ny
    subplot(ny ,1,idx)
    plot(t,y(:,idx),'Linewidth',1);  
    if idx == 1
        title('WECS outputs');
    end
    axis tight;
    grid on;
    ylabel(sysOut.OutputName{idx});
end
xlabel('time [s]');
cl = colormap('lines');

figure(2)
idxEnd = size(x,2)-1;
for idx = 1:idxEnd
    subplot(idxEnd,1,idx)
    plot(t,x(:,idx),'-','Color',cl(idx,:));
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

plot(t, d,'k--','Linewidth',1);
hold on;
plot(t,x(:,idx),'-','Color',cl(1,:));
plot(t,xest(:,idx),'--','Color',cl(2,:));
legend('wind','wind state','wind est');
ylabel(' Wind [m/s]')
xlabel('time [s]');

