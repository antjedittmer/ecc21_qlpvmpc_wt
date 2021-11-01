temp1=fft(OutDataWind(:,2));
Fs = 1/DT;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(temp1);             % Length of signal
t = (0:L-1)*T;  
P2 = abs(temp1/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);% Time vector
f = Fs*(0:(L/(2)))/L;
figure;plot(f(1:100),P1(1:100)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
temp2=fft(OutDataWind(:,1));
Preal2 = abs(temp2/L);
Preal1 = Preal2(1:L/2+1);
Preal1(2:end-1) = 2*Preal1(2:end-1);
figure;plot(f(1:100),P1(1:100),f(1:100),Preal1(1:100))
ylim