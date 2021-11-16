clear; clc; bdclose all; close all; 
addpath(fullfile(fileparts(pwd),'datain'));
mdlName = 'test_SimulinkMdl2_qLPVMPCbeta_UIO_testOL'; %
% load('workSpace.mat')
% load('OutDataSweep.mat');
% 
% param;
% initModel5MWNREL;
% initModel5MWNRELBianchi;

exUIO_WECS;
close all;

legStates = {'Torque_g_r [kNm]','\beta_r [rad]','omega_r [rad/s]', 'omega_g [rad/s]', 'Torque_g [kNm]','\beta [rad]',};
ysInputname = {'V [m/s]', 'Torque_g [kNm]', '\beta_d [rad]',};


indexVec = 1:22;
lenRun = length(indexVec);
yoff = nan(5,lenRun);
xoff = nan(length(legStates),lenRun);
xoffEst = nan(5,lenRun);
uoff = nan(3,lenRun);

Vbar = 4;
for index = indexVec
    Vbar = Lin_points.V(index);
    omegabar = Lin_points.RSpeed(index);
    Lambda_bar = wecs.rb*omegabar./ Vbar;
    
%     beta_bar = Lin_points.Pitch(index)/pi*180;
%     Cqb = interp2(Rotor_Lambda,Rotor_Pitch,Cq,Lambda_bar,beta_bar,'cubic');
%     
    %Torque_g_bar = 0.5*rho*pi*wecs.rb^3*Cqb*Vbar^2/G;
    simOut = sim(mdlName);
    
    xoff(:,index) = OutDataTest1(end,:);
    xoffEst(:,index) = OutDataTestLin(end,:);
%     xoff(:,index) = [[q(end,1:2), q(end,3)-q(end,4)], qdot(end,:), beta(end)]';
%     uoff(:,index)= [Vbar, beta_bar/180*pi,Torque_g_bar]';
end

n = 10;
figure(1+n);
for idx = 1: 3
    subplot(3,1,idx);
    plot(indexVec +3,uoff(idx,:))
    axis tight; grid on;
    ylabel(sysInputname{idx});
end


figure(2+n);
for idx = 1: 4
    
    subplot(4,1,idx);
    plot(indexVec +3,xoff(idx,:))
     axis tight; grid on;
    ylabel(legStates{idx})
end

figure(3+n);
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
