clear; clc; bdclose all; close all; 
addpath(fullfile(fileparts(pwd),'datain'));
mdlName = 'testBianchi.slx'; %
load('workSpace.mat')
load('OutDataSweep.mat');

param;
initModel5MWNREL;
initModel5MWNRELBianchi;

legStates = {'yt' 'zeta' 'thetha_s' 'yt_dot' 'zeta_dot' 'omega_r' 'omega_g' 'beta'};
 sysInputname = {'V [m/s]', '\beta_d [rad]', 'Torque_g [kNm]'};


indexVec = 1:22;
lenRun = length(indexVec);
yoff = nan(2,lenRun);
xoff = nan(length(legStates),lenRun);
uoff = nan(3,lenRun);

for index = indexVec
    Vbar = Lin_points.V(index);
    omegabar = Lin_points.RSpeed(index);
    Lambda_bar = omegabar./ Vbar;
    
    beta_bar = Lin_points.Pitch(index)/pi*180;
    Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
    
    Torque_g_bar = 0.5*rho*pi*rb^3*Cqb*Vbar^2/G;
    sim(mdlName);
    
    yoff(:,index) = [qdot(end,end)*G, Tg(end)]';
    xoff(:,index) = [[q(end,1:2), q(end,3)-q(end,4)], qdot(end,:), beta(end)]';
    uoff(:,index)= [Vbar, beta_bar/180*pi,Torque_g_bar]';
end

figure(1);
for idx = 1: 3
    subplot(3,1,idx);
    plot(indexVec +3,uoff(idx,:))
    axis tight; grid on;
    ylabel(sysInputname{idx});
end


figure(2);
for idx = 1: 4
    
    subplot(4,1,idx);
    plot(indexVec +3,xoff(idx,:))
     axis tight; grid on;
    ylabel(legStates{idx})
end

figure(3);
for idx = 5: 8
    subplot(4,1,idx-4);
    plot(indexVec +3,xoff(idx,:))
    axis tight; grid on;
    ylabel(legStates{idx})
end


save('linsys_LPV',  'xoff', 'yoff', 'uoff'); % , 'omega_sched', 'U_sched')



% %% Original code
% nULin = length(Vbar);
% ULinspace = lVbar;
% 
% params(1).Name = 'ULin';
% params(1).Value = ULingrid;
% %     params(2).Name = 'tauLin';
% %     params(2).Value = tauLingrid;
% 
% %%
% op = findop(mdl, 250, params);
% 
% %%
% io = getlinio(mdl);
% tmp_io = io(4,1);
% io(4,1) = io(5,1);
% io(5,1) = tmp_io;
% 
% %%
% [linsys, ~, info] = linearize(mdl, op, io, params, ...
%     linearizeOptions('StoreOffsets', true));
% 
% %%
% offsets = getOffsetsForLPV(info);
% yoff = offsets.y;
% xoff = offsets.x;
% uoff = offsets.u;
% 
% omega_sched = reshape(squeeze(yoff(1,:,:)), [nULin ntauLin]);
% U_sched = ULinspace;
% 
% linsys.SamplingGrid = struct('omega_r',omega_sched);
% 
% linsys.InputName = {'Wind', 'Torque'};
% linsys.OutputName = {'omega_r', 'A_y', 'P_g'};
% linsys.StateName = {'q1', 'q2', 'q3', 'q4', 'omega_r'};
% 
% save('linsys_LPV', 'linsys', 'xoff', 'yoff', 'uoff', 'omega_sched', 'U_sched')
