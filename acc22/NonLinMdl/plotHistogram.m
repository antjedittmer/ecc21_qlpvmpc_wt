function plotHistogram(OutTableTest2Plot, OutTableMPCPlot,OutTableMPC2Plot,PGRef)

% Set path to figure directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);
figDir = fullfile(mainDir,'figDirUIO');

subplot(1,2,1);
histogram(OutTableTest2Plot.NcIMUTAxs,20); hold on;
histogram(OutTableMPCPlot.NcIMUTAxs,20); histogram(OutTableMPC2Plot.NcIMUTAxs,20); 
hold off;
axis tight; grid on;
ylabel('Occurence [-]');
xlabel( 'Twr_{FA} y_t [m/s^2]')
subplot(1,2,2);
% histogram(abs(OutTableTest2Plot.GenPwr - PGRef),20); hold on;
% histogram(abs(OutTableMPCPlot.GenPwr - PGRef),20); histogram(abs(OutTableMPC2Plot.GenPwr - PGRef),20); 
% xlabel('|P_{g,ref}-P_g| [kW]')
histogram(OutTableTest2Plot.GenPwr - PGRef,20); hold on;
histogram(OutTableMPCPlot.GenPwr - PGRef,20); histogram(OutTableMPC2Plot.GenPwr - PGRef,20); 
xlabel('P_{g,ref}-P_g [kW]')
legend('P-PI','MPC V_{est}','MPC V_{real}','Position',[0.3285 0.7738 0.1792 0.1341])
hold off;
axis tight; grid on;
% save to png/eps
print(gcf,fullfile(figDir,'histTwrPabs'), '-dpng');
print(gcf,fullfile(figDir,'histTwrPabs'), '-depsc');
