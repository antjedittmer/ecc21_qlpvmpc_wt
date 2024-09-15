varnames = {'wind', 'RotSpd', 'GenSpd', 'Pitch', 'Tg_meas',  ...
    'ytdodot','xdotdot','Theta_s', 'omega_gr'};

OutTableTest1 = array2table(OutDataTest,'VariableNames',varnames);
OutTableTest1.GenPwr = OutTableTest1.GenSpd  .* OutTableTest1.Tg_meas * 0.944/1000; 

OutTableTestFAST = array2table(OutDataTestFAST,'VariableNames',varnames);
OutTableTestFAST.GenPwr = OutTableTestFAST.GenSpd  .* OutTableTestFAST.Tg_meas * 0.944/1000; 

time = 0 : DT : (length(OutDataTest)*DT - DT);

figure(1);
ax1 = nan(9,1);
for idx = 1:9
    ax1(idx) = subplot(9,1,idx);
    plot(time,OutTableTest1{:,idx},time, OutTableTestFAST{:,idx},'--' );
    ylabel(varnames{idx})
    axis tight; grid on;
end
xlabel('time [s]')
legend('Simulink','FAST');
linkaxes(ax1,'x')

% 
% posDefault = get(gcf, 'position');
% set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.5]);


figure(2);
plot(time, OutTableTest1.GenPwr, time, OutTableTestFAST.GenPwr,'--');