% rga.m      
%
% This function calculates the relative gain array for a transfer function
% matrix of a 2X2 LTI system within given boundaries
% 
%
% Syntax:  lambda = rga(G, lowerb, upperb)
%
% Input parameters:
%    G         - Transfer function matrix of a 2X2 LTI system
%    lowerb    - 10^lowerb 1/s is the lower frequency boundary
%    upperb    - 10^upperb 1/s is the upper frequency boundary
%
%
% Output parameters:
%    lambda    - lambda(:,:,k) is the rga of G at the frequency point k
%                
%
% Author:     Birgit Menzel

function lambda = plot_rga(G,lowerb,upperb,titleStr)

if nargin < 4
    titleStr = '';
end

SamplePoints = (upperb-lowerb)*1000;
s= 1i *logspace(lowerb,upperb,SamplePoints);
for k=1:1:SamplePoints
    Gsample(:,:,k)=evalfr(G,s(1,k));
    lambda11(1,k)= Gsample(1,1,k) * Gsample(2,2,k)/(Gsample(1,1,k)*Gsample(2,2,k)-Gsample(1,2,k)*Gsample(2,1,k));
    lambda12(1,k)= -Gsample(1,2,k) * Gsample(2,1,k)/(Gsample(1,1,k)*Gsample(2,2,k)-Gsample(1,2,k)*Gsample(2,1,k));
    lambda21(1,k)= lambda12(1,k);
    lambda22(1,k)= lambda11(1,k);
    lambda(:,:,k)= [lambda11(1,k) lambda12(1,k); lambda21(1,k) lambda22(1,k)];
end
    
f=abs(s);
lambda11(1,:)=lambda(1,1,:);
lambda12(1,:)=lambda(1,2,:);
lambda21(1,:)=lambda(2,1,:);
lambda22(1,:)=lambda(2,2,:);

% g11(1,:)=abs(Gsample(1,1,:));
% g12(1,:)=abs(Gsample(1,2,:));
% g21(1,:)=abs(Gsample(2,1,:));
% g22(1,:)=abs(Gsample(2,2,:));

% plot of the real and imaginary parts of the rga over frequency
figure;
style=12;
set(0,'defaultaxesfontsize',style,'defaulttextfontsize',style,...
    'defaultlineLinewidth',1.5)
title('Relative Gain Array of the Transfer Function in Operating Point No 5');
subplot(2,2,1);
semilogx(f,real(lambda11), 'k-');
hold on;
semilogx(f,imag(lambda11), 'k:');
xlabel('Frequency [rad]');
ylabel('RGA_{11}');
title(titleStr)
grid on; axis tight;
hold off;

subplot(2,2,2);
semilogx(f,real(lambda12), 'k-');
hold on;
semilogx(f,imag(lambda12), 'k:');
xlabel('Frequency [rad]');
ylabel('RGA_{12}');
grid on; axis tight;
hold off;

subplot(2,2,3);
semilogx(f,real(lambda21), 'k-');
hold on;
semilogx(f,imag(lambda21), 'k:');
xlabel('Frequency [rad]');
ylabel('RGA_{21}');
grid on; axis tight;
legend('Real','Imag','Location','NorthWest');
hold off;

subplot(2,2,4);
semilogx(f,real(lambda22),'k-');
hold on;
semilogx(f,imag(lambda22),'k:');
xlabel('Frequency [rad]');
ylabel('RGA_{22}');
grid on; axis tight;
hold off;