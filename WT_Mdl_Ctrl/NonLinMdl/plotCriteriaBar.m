function plotCriteriaBar(K_matrix,strK,summaryTickLabel)
figure;
cl1 = colormap('lines');
cl = [0.5,0.5,0.5;
    cl1(4,:);
    [0.8,0,0];
    cl1(1,:);
    cl1(5,:)];

b = bar(K_matrix);
legend('C_{P}','C_{\sigma(P)}','C_{Actuator}','C_{D,Blades}','C_{D,Tower}','Location','NorthEast');

for k = 1:5
    b(k).FaceColor = cl(k,:);
end

axis tight; grid on;
% titleCellNew = {['Power tracking C_{P}', '{\color[rgb]{',num2str(cl(1,:)),'}',strK.power(1:end-2),'}'],...
%     ['Power variance C_{\sigma(P)}: ',  '{\color[rgb]{',num2str(cl(2,:)),'}',strK.StdPowSpeed(1:end-2),'}'], ...
%     ['Actuator power C_{Actuator}: ',  '{\color[rgb]{',num2str(cl(3,:)),'}',strK.ActPwr(1:end-2),'}'], ...
%     ['Blade damage C_{D,Blades}: ',  '{\color[rgb]{',num2str(cl(4,:)),'}',strK.DamageBlades(1:end-2),'}'],...
%     ['Tower damage C_{D,Tower}: ',  '{\color[rgb]{',num2str(cl(5,:)),'}',strK.DamageTower(1:end-2),'}']};

%title(titleCellNew)

ax = gca;
ax.XTickLabel = summaryTickLabel;