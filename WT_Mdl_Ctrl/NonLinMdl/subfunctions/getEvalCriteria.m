function [K,W,ref] = getEvalCriteria(aTable,ref)
% calculate evaluation criteria for controllers

% epected variablenames: GenPwr,GenSpeed,TwrBsM,BlPitchC
if nargin <2 || isempty(ref)
        ref = [];
end

% Get quantitative criteria for controller evaluation
idxTime = aTable.Time > 15; % 10 : height(outTable1); %;  %
varNames = aTable.Properties.VariableNames;


%% Get evaluation criteria
%Get mean generator power
meanGenPwr = mean(aTable.GenPwr(idxTime)/1000);

%Get std generator power
stdGenPwr = std(aTable.GenPwr(idxTime));

%Get std generator speed power
%stdGenSpeed = std(aTable.GenSpeed(idxTime));

%C_damageBlades calculated via Rainflow count w flap. blade RootMyb signals
% Antje ToDo Using flapping moments only: Also edgewise(x) and pitch (z)?
% Antje ToDo: Transform from rotating into fixed reference system?
idxRootMyb = contains(varNames,'RootMyb') | contains(varNames,'RootMFlp');
m = 10;
damageBlades = calculateDamageBlades(aTable{idxTime,idxRootMyb},m);


%C_damageTower calculated via Rainflow count with tower base moments 
idxTowerM = contains(varNames,'TwrBsM');
mT = 4;
damageTower = calculateDamageBlades(aTable{idxTime,idxTowerM},mT);


% C_actuator is the magnitude of the power consumed by the actuators
idxBlPitchC = contains(varNames,'BlPitch');
 %flapwise moment (i.e., the moment caused by flapwise forces) at the blade root
% AD ToDo: Should this be in blade/rotating frame or fixed? Definition P_{act,ref} = P_{act,0} 
meanActPwr = mean(mean(abs(diff(aTable{idxTime,idxBlPitchC}).*aTable{idxTime,idxRootMyb}(1:end-1,:))));


% Opt Criteria
if isempty(ref)
    ref.Pw = 5;
    %ref.StdGenSpeed = stdGenSpeed;% max(stdGenSpeed);
    ref.StdGenPwr = stdGenPwr;% max(stdGenPwr);
    ref.MeanActPwr = meanActPwr;% max(meanActPwr); %Definition P_{act,ref} = P_{act,0} ?
    ref.DamageBlades = damageBlades;%max(damageBlades);
    ref.DamageTower = damageTower;%max(damageTower);
    K = NaN;
    W = NaN;
    
else
    K.GenPwr  = abs(meanGenPwr/ref.Pw - 1);
    %K.stdGenSpeed = stdGenSpeed/ref.StdGenSpeed;
    K.stdGenPwr = stdGenPwr/ref.StdGenPwr;
    K.ActPwr = meanActPwr/ref.MeanActPwr;
    K.damageBlades = damageBlades/ref.DamageBlades;
    K.damageTower = damageTower/ref.DamageTower;
    W = NaN;
    
end
