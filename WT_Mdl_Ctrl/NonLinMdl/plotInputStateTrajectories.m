function plotInputStateTrajectories(Uvec,xest,time, titleStr1, titleStr2,nrPlot1)

N = size(Uvec,2)/2;
Nsim = length(time);
DtPlot = round(Nsim/5);

xCell = {'omega_g [rad]','T_g','ThetaDot_s [rad/s]','omegaDot_r [rad/s^2]',...
    'omegaDot_gr [rad/s^2]','TDot_g'}; %'V_dot'
uCell = {'betaDot [rad/s]','TgDot [Nm/s]'}; %'V_dot'
legCell = regexp(deblank(sprintf('time:%2.3fs ', time(1:DtPlot:Nsim))),' ','split');



figure(nrPlot1);
for idx = 1:2
    for idx1 = 1:DtPlot:Nsim
        subplot(2,1,idx);
        plot(1:N,Uvec(idx1,idx:2:end)); hold on; axis tight; grid on
    end
    legend(legCell)
    if idx == 1
        title(titleStr1)
    end
    axis tight; grid on;
    ylabel(uCell{idx});
    hold off;
end


figure(nrPlot1+1);
for idx = 1:6
    subplot(6,1,idx);
    plot(xest(:,idx))
    axis tight; grid on;
    legend(xCell{idx});
    if idx == 1
        title(titleStr2)
    end
end