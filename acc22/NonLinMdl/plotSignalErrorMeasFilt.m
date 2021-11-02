function plotSignalErrorMeasFilt(OutDataRotorSpeed,DT,ylabStr)
% plotSignalErrorMeasFilt plots the filtered and measured signal erro

%% Assign default inputs
if nargin <1 || isempty(OutDataRotorSpeed)
    load(fullfile('dataOut','OutDataRotorSpeedPI'),'OutDataRotorSpeed');
end

if nargin <2 || isempty(DT)
    DT = 0.008;
end
if nargin <3 || isempty(ylabStr)
    ylabStr = 'Rotor Speed [rad/s]';
end

%% Create signals with names for plots

% Time and signal values
timeVecAll = (0:DT:(length(OutDataRotorSpeed)-1)*DT)'; % time

idx = timeVecAll > 15;
timeVecData = timeVecAll(idx);

rotSpd = OutDataRotorSpeed(idx,1);
filtRotSpd = OutDataRotorSpeed(idx,3);
measRotSpd = OutDataRotorSpeed(idx,2);

% error signals
errFiltRotSpd = rotSpd - filtRotSpd;
errMeasRotSpd = rotSpd-  measRotSpd;

% string with std and VAF for plots
str1 = getStdVAFstring(errFiltRotSpd,rotSpd);
str2 = getStdVAFstring(errMeasRotSpd,rotSpd);

cl = colormap('lines');
figure;
subplot(2,2,1)
plot(timeVecData,rotSpd,'-',timeVecData,filtRotSpd,'--');
legend('True','Filt'); 
axis tight; grid on;
ylabel(ylabStr);

subplot(2,2,2)
plot(timeVecData,errFiltRotSpd,'Color',cl(2,:));
axis tight; grid on;
legend(['True - Filt',str1]); 
axis tight; grid on;

subplot(2,2,3)
plot(timeVecData,measRotSpd,'--','Color',cl(3,:,:));
hold on;
plot(timeVecData,rotSpd,'-','Color',cl(1,:,:));
legend('Meas','True'); 
axis tight; grid on;
ylabel(ylabStr);

subplot(2,2,4)
plot(timeVecData,errMeasRotSpd,'Color',cl(3,:));
axis tight; grid on;
legend(['True - meas',str2]); 
axis tight; grid on;

function str1 = getStdVAFstring(errFiltRotSpd,rotSpd)
stdFilt = std(errFiltRotSpd);
vafFilt = 1 - var(errFiltRotSpd)/var(rotSpd);
str1 = sprintf(': std: %2.2e, VAF: %d%%',stdFilt, round(vafFilt*100));

