% clc;
% close all;
tau_0=0;
computingInstant=1;
paths=4;
% BETA=rand(1,paths);
% TAUJ=(rand(1,paths)*2-1);
rxsignal=[];
BETA=[0.2691 0.4228 0.5479 0.9247];
TAUJ=[0.0042 0.0098 0.0030 0.0070];
z=1;
t1=1;

df = 1/1000;
freqRange = 0:df:1;
transferFunction=zeros(length(freqRange), 1);

i = 1;
for f=freqRange
    temp1=0;
    for p=1:1:paths
        temp1=temp1+BETA(p)*exp(-1i*2*pi*f*tau_0)*exp(-1i*2*pi*f*TAUJ(p)*computingInstant);
    end
    transferFunction(i) = temp1;
    i = i + 1;
end

figure;
plot(freqRange*1000,abs(transferFunction));
title(sprintf('Time Varying Transfer Function at the time instan t=%dus', t0))








%% backup

clc
clear all
close all

Tau0 =0;
f=1;
nop =4;
% BETA = rand(1,nop);
Tauj = [0.62 1.84 0.86 0.37];
BETA = [0.23 0.17 0.23 0.44];
rxsignal =[];
tvtf = [];
t = 0:(1/100):99.99;
txsignal = cos(2*pi*f*t);
z=1;
for t =0:(1/100):99.99
    temp =0;
    temp1 =0;
    for p=1:1:nop
        temp = temp+BETA(p)*exp(2*pi*f*t)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(p)*t);
        temp1 = temp1+BETA(p)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(p)*t);
    end
    rxsignal = [rxsignal temp];
    tvtf = [tvtf temp1];
  
end

figure(1)
subplot(2,2,1)
plot(txsignal)
axis([1 1000 -2 2]);
title('Transmitted signal');

subplot(2,2,2)
fre = (0:1:length(rxsignal)-1)/100;
plot(real(rxsignal),'r')
axis([1 1000 -2 2]);
title('Real part of received signal');

subplot(2,2,3)
plot(fre, abs(fft(txsignal)));
axis([0 2 0 1000]);
title('Spectrum of transfer signal');

subplot(2,2,4)
plot(fre, abs(fft(real(rxsignal))));
title('real part of corresponding spectrum of transmitted signal');

figure(2)
subplot(2,1,1)
plot(abs(tvtf));
axis([0 1000 0 2]);
title('Time varying transfer function value magnitude');

subplot(2,1,2)
plot(phase(tvtf));
axis([0 1000 -25 0]);
title('Time varying transfer function value phase');

%% raghav (fast fading mayebe)
clc;
clear all;
close all;

t=0:2:100;
fo = 1;
%J=No of Paths
J=4;
B = [0.23 0.17 0.23 0.44];
T = [0.62 1.84 0.86 0.37];

D=T.*(-fo);
fspread = D + fo;

transmitted = cos(2*pi*fo*t);
y=zeros(1,length(t));
H = zeros(1,length(t));
for k=1:length(t)
    for j=1:J
        y(k)=y(k)+B(j)*exp(1i*2*pi*fspread(j)*t(k))*exp(-1i*2*pi*fo);
        H(k)=H(k)+B(j)*exp(1i*2*pi*D(j)*t(k))*exp(-1i*2*pi*fo)*exp(i*2*pi*fo*t(k));
    end
end

figure(1);
subplot(2,2,1);
plot(transmitted);
title('transmitted signal');
xlabel('time');
ylabel('Amplitude');

subplot(2,2,3);
plot(abs(fft(transmitted)));
title('Spectrum transmitted signal');
xlabel('frequency');
ylabel('Amplitude');

subplot(2,2,2);
plot(t,real(y));
title('received signal');
xlabel('time');
ylabel('Amplitude');

subplot(2,2,4);
plot(abs(fft(y)));
title('Spectrum of received signal');
xlabel('time');
ylabel('Amplitude');

figure(2);
subplot(2,1,1);
plot(abs(H));
title('Magnitude of Fast-Fading Transfer function time varying Doppler Spread signal');
xlabel('time');
ylabel('Amplitude');

subplot(2,1,2);
plot(-phase(H));
title('Phase of Fast-Fading Transfer function time varying Doppler Spread signal');
xlabel('time');
ylabel('Angle');