function calculateEvalCriteria(cellOfTable,summaryTickLabel)
% calculate evaluation criteria for controllers

% epected variablenames: GenPwr,GenSpeed,TwrBsM,BlPitchC
outTable1 = cellOfTable{1};
lenComparison = length(cellOfTable);
% Get quantitative criteria for controller evaluation
idxTime = outTable1.Time > 15; % 10 : height(outTable1); %;  %
varNames = outTable1.Properties.VariableNames;

isFAST = any(contains(varNames,'GenSpeed'));

if nargin <2
    summaryTickLabel = {'CPC', 'IPC_{fminsearch}','IPC_{SimAnn}'};
    %summaryTickLabel = {'NREL_{FAST}', 'NREL_{OC3}','DLR','TUDelft'};
end

%% Get evaluation criteria

%Get mean generator power
meanGenPwr = nan( 1,lenComparison);
for idx = 1: lenComparison
    meanGenPwr(idx) = mean(cellOfTable{idx}.GenPwr(idxTime)/1000);
end

%Get std generator power
stdGenPwr = nan( 1,lenComparison);
for idx = 1: lenComparison
    stdGenPwr(idx) = std(cellOfTable{idx}.GenPwr(idxTime));
end

%Get std generator speed power

stdGenSpeed = nan( 1,lenComparison);
if isFAST
for idx = 1: lenComparison
    stdGenSpeed(idx) = std(cellOfTable{idx}.GenSpeed(idxTime));
end
end

%C_damageBlades calculated via Rainflow count w flap. blade RootMyb signals
% Antje ToDo Using flapping moments only: Also edgewise(x) and pitch (z)?
% Antje ToDo: Transform from rotating into fixed reference system?
 if isFAST
idxRootMyb = contains(varNames,'RootMyb') | contains(varNames,'RootMFlp');
if any(idxRootMyb)
    m = 10;
    damageBlades = nan( 1,lenComparison);
    for idx = 1: lenComparison
        damageBlades(idx) = calculateDamageBlades(cellOfTable{idx}(idxTime,:),m,idxRootMyb);
    end
end

 else
    damageBlades = nan( 1,lenComparison);
    for idx = 1: lenComparison
        damageBlades(idx) = std(cellOfTable{idx}.zeta_dotdot(idxTime,:));
    end
 end

% C_damageTower calculated via Rainflow count with tower base moments
 if isFAST
idxTowerM = contains(varNames,'TwrBsM');
if any(idxTowerM) % Antje ToDo What to use for offshore/floating wind turbines?
    m = 4;
    damageTower = nan( 1,lenComparison);
    for idx = 1: lenComparison
        damageTower(idx) = calculateDamageBlades(cellOfTable{idx}(idxTime,:),m,idxTowerM);
    end
else
    damageTower = nan(size(damageBlades));
end
 else
    damageTower = nan( 1,lenComparison);
    for idx = 1: lenComparison
        damageTower(idx) = std(cellOfTable{idx}.NcIMUTAxs(idxTime,:));
    end
 end

% C_actuator is the magnitude of the power consumed by the actuators
idxBlPitchC = contains(varNames,'BlPitch');
 %flapwise moment (i.e., the moment caused by flapwise forces) at the blade root
% AD ToDo: Should this be in blade/rotating frame or fixed? Definition P_{act,ref} = P_{act,0} 
meanActPwr = nan( 1,lenComparison);
 if isFAST
for idx = 1: lenComparison
      meanActPwr(idx) = ...
          mean(mean(abs(diff(cellOfTable{idx}{idxTime,idxBlPitchC}).*cellOfTable{idx}{idxTime,idxRootMyb}(1:end-1,:))));
end
 else
     for idx = 1: lenComparison
      meanActPwr(idx) = ...
          mean(abs(diff(cellOfTable{idx}{idxTime,idxBlPitchC})));
     end
end

% Opt Criteria
refPw = abs(meanGenPwr(1)/5 - 1);
refStdGenSpeed = stdGenSpeed(1); % max(stdGenSpeed);
refStdGenPwr = stdGenPwr(1); % max(stdGenPwr);
refMeanActPwr = meanActPwr(1); % max(meanActPwr); %Definition P_{act,ref} = P_{act,0} ?
refDamageBlades = damageBlades(1); %max(damageBlades);
refDamageTower = damageTower(1); %max(damageTower);

K_GenPwr  = abs(meanGenPwr/refPw - 1);
K_stdGenSpeed = stdGenSpeed/refStdGenSpeed;
K_stdGenPwr = stdGenPwr/refStdGenPwr;
K_ActPwr = meanActPwr/refMeanActPwr;
K_damageBlades = damageBlades/refDamageBlades;
K_damageTower = damageTower/refDamageTower;

cl = colormap('lines');

strMeanGenPwr = sprintf('%2.3f, ', meanGenPwr);
strKpowerabs = sprintf('%2.3f, ', K_GenPwr);
%strStdGenSpeed = sprintf('%2.3f, ', stdGenSpeed);

strKpower = sprintf('%2.3f, ', K_GenPwr/max(K_GenPwr));
strKStdGenSpeed = sprintf('%2.3f, ',K_stdGenSpeed);
strKStdPowSpeed = sprintf('%2.3f, ',K_stdGenPwr);
strKActPwr = sprintf('%2.3f, ',K_ActPwr);
strKDamageBlades = sprintf('%2.3f, ',K_damageBlades);
strKDamageTower = sprintf('%2.3f, ',K_damageTower);

figure;


bar([K_GenPwr/max(K_GenPwr);K_stdGenSpeed]');
axis tight; grid on;
titleCell = {'Comparison PI Torque/CPC Controller',...
    ['Mean GenPwr [MW]: ',strMeanGenPwr(1:end-2)],...
    ['C_{Power,abs} = |mean(P_{ctrl})/mean(P_{max}) - 1|: ', strKpowerabs(1:end-2)],...
    ['C_{Power} = C_{Power,abs}/C_{Power,amax}: ', '{\color[rgb]{',num2str(cl(1,:)),'}',strKpower(1:end-2),'}'],...
    ['C_{StdGenSpeed} = |STD_{ctrl}/STD_{max}|: ',  '{\color[rgb]{',num2str(cl(2,:)),'}',strKStdGenSpeed(1:end-2),'}']};
title(titleCell)
xlabel('Controller');
ylabel('Minimization criteria [-]');
legend('C_{Power}','C_{StdGenSpd}','Location','SouthEast')
ax = gca;
ax.XTickLabel = summaryTickLabel;

figure;
bar([K_GenPwr/max(K_GenPwr); K_stdGenPwr; K_ActPwr; K_damageBlades; K_damageTower]');
axis tight; grid on;
titleCellNew = {titleCell{4},...
    ['C_{StdPwr} = |STD_{ctrl}/STD_{ref}|: ',  '{\color[rgb]{',num2str(cl(2,:)),'}',strKStdPowSpeed(1:end-2),'}'], ...
    ['C_{ActPwr} = |ActPwr_{ctrl}/ActPwr_{ref}|: ',  '{\color[rgb]{',num2str(cl(3,:)),'}',strKActPwr(1:end-2),'}'], ...
    ['C_{D,Blades} = |DBlades_{ctrl}DBlades_{ref}|: ',  '{\color[rgb]{',num2str(cl(4,:)),'}',strKDamageBlades(1:end-2),'}'],...
    ['C_{D,Tower} = |DBlades_{ctrl}DBlades_{ref}|: ',  '{\color[rgb]{',num2str(cl(5,:)),'}',strKDamageTower(1:end-2),'}']};

title(titleCellNew)
ax = gca;
ax.XTickLabel = summaryTickLabel;