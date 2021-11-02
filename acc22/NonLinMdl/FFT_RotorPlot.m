clc; close all; clear;

matFileRotorSpeed = fullfile('dataOut','OutDataRotorSpeedMPC.mat');
DT = 0.008;

load(matFileRotorSpeed,'OutDataRotorSpeed');
S = OutDataRotorSpeed(:,1);
S2 = OutDataRotorSpeed(:,3);
S3 = OutDataRotorSpeed(:,2);

figure;
timeVecAll = (0:DT:(length(OutDataRotorSpeed)-1)*DT)';
subplot(2,1,1)
plot(timeVecAll,S3,timeVecAll,S2,timeVecAll,S);

% rotSpd = OutDataRotorSpeed(idx,1);
% filtRotSpd = OutDataRotorSpeed(idx,3);
% measRotSpd = OutDataRotorSpeed(idx,2);

Fs = 1/DT;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(S); % Length of signal
t = (0:L-1)*T;        % Time vector

Ldot5 = round(L/2);
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:Ldot5+1);
P1(2:end-1) = 2*P1(2:end-1);

Y2 = fft(S2);
P2_2 = abs(Y2/L);
P1_2 = P2_2(1:Ldot5+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

Y3 = fft(S3);
P2_3 = abs(Y3/L);
P1_3 = P2_3(1:Ldot5+1);
P1_3(2:end-1) = 2*P1_3(2:end-1);

f = Fs*(0:Ldot5)/L;
idxF = f<1;
idxEnd = sum(idxF);


P1cumSum = cumsum(P1(2:idxEnd));
idxF = min(find(P1cumSum > 0.95 * sum(P1(2:idxEnd))));

minF2 = 3/(2*pi);
idxF2 = find(min(abs(f - minF2)) == abs(f-minF2));
percentP1 = P1cumSum(idxF2)/sum(P1(2:idxEnd));

figure;
plot(f(2:idxEnd),P1(2:idxEnd),...
    f(2:idxEnd), P1_2(2:idxEnd),...
    f(2:idxEnd), P1_3(2:idxEnd),...
    f(idxF),P1(idxF),'g*',...
    f(idxF2),P1(idxF2),'k*'); 
legend('True','Filt','Meas','f = 0.84 (99%)','f = 0.48 (94%)');
%title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
