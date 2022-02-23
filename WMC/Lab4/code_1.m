%% fast fading

clc;
clear all;
close all;
 
Tau0 =0;
f=1;
nop =4;
Tauj = [0.62 1.84 0.86 0.37];
BETA = [0.23 0.17 0.23 0.44];
fshift=[];
for j=1:1:nop
    fshift(j)=abs(-f+Tauj(j));
end
rxsignal =[];
tvtf = [];
t = 0:(1/100):100;
txsignal = cos(2*pi*f*t);
for t =0:(1/100):100
    rx =0;
    tf =0;
    for p=1:1:nop
        rx = rx+ BETA(p)*cos(2*pi*fshift(p)*t);
        tf = tf+BETA(p)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(p)*t);
    end
    rxsignal = [rxsignal rx];
    tvtf = [tvtf tf];
end
 
figure(1)
subplot(2,2,1)
plot(txsignal)
axis([1 1000 -2 2]);
title('Transmitted signal U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('amplitude');
 
subplot(2,2,2)
fre = (0:1:length(rxsignal)-1)/100;
plot(real(rxsignal),'r')
axis([1 1000 -2 2]);
title('Received signal U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('amplitude');
 
subplot(2,2,3)
plot(fre, abs(fft(txsignal)));
axis([0 2 0 1000]);
title('Spectrum of transfer signal U19EC046');
xlabel('Frequency in MHz');
ylabel('amplitude');
 
subplot(2,2,4)
plot(fre, abs(fft(real(rxsignal))));
axis([0 2 0 1000]);
title('Corresponding spectrum of received signal U19EC046');
xlabel('Frequency in MHz');
ylabel('amplitude');
 
figure(2)
subplot(2,1,1)
plot(abs(tvtf));
axis([0 1000 0 2]);
title('Time varying transfer function value magnitude U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('amplitude');
 
subplot(2,1,2)
plot(phase(tvtf));
axis([0 1000 -25 0]);
title('Time varying transfer function value phase U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('angle');

%% slow fading

clc;
clear all;
close all;
 
Tau0 =0;
f=1;
nop =4;
 
Tauj = [0.0042 0.0098 0.0030 0.007];
BETA = [0.2691 0.4228 0.5479 0.9427];
fshift=[];
for j=1:1:nop
    fshift(j)=abs(-f+Tauj(j));
end
rxsignal =[];
tvtf = [];
t = 0:(1/100):100;
txsignal = cos(2*pi*f*t);
for t =0:(1/100):100
    temp =0;
    temp1 =0;
    for p=1:1:nop
        temp = temp+ BETA(p)*cos(2*pi*fshift(p)*t);
        temp1 = temp1+BETA(p)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(p)*t);
    end
    rxsignal = [rxsignal temp];
    tvtf = [tvtf temp1];
end
 
figure(1)
subplot(2,2,1)
plot(txsignal)
axis([1 10000 -2 2]);
title('Transmitted signal U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('amplitude');
 
subplot(2,2,2)
fre = (0:1:length(rxsignal)-1)/100;
plot(real(rxsignal),'r')
axis([1 10000 -2 2]);
title('Real part of received signal U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('amplitude');
 
subplot(2,2,3)
plot(fre, abs(fft(txsignal)));
axis([0 2 0 1000]);
title('Spectrum of transfer signal U19EC046');
xlabel('Frequency in MHz');
ylabel('amplitude');
 
subplot(2,2,4)
plot(fre, abs(fft(real(rxsignal))));
axis([0 2 0 1000]);
title('Real part of corresponding spectrum of received signal U19EC046');
xlabel('Frequency in MHz');
ylabel('amplitude');
 
figure(2)
subplot(2,1,1)
plot(abs(tvtf));
axis([0 10000 0 2.3]);
title('Time varying transfer function value magnitude U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('amplitude');
 
subplot(2,1,2)
plot(phase(tvtf));
axis([0 10000 -4 0]);
title('Time varying transfer function value phase U19EC046');
xlabel('Samples in time domain with Ts = 1/100 micro secs');
ylabel('angle');

%% delay spread

clc;
clear all;
close all;
 
TAU0=0;
t0=1;
nop=4;
BETA=rand(1,nop);
TAUJ=(rand(1,nop)*2-1);
rxsignal=[];
BETA= [0.9575 0.9649 0.1576 0.9706];
TAUJ= [0.9143 -0.0292 0.6006 -0.7162];
tv_tf_comp_at_t0=[];
z=1;
t1=1;
 
for f=0:(1/1000):0.999
    temp=0;
    temp1=0;
    for p=1:1:nop
        temp1=temp1+BETA(p)*exp(-1i*2*pi*f*TAU0)*exp(-1i*2*pi*f*TAUJ(p)*t0);
    end
    tv_tf_comp_at_t0=[tv_tf_comp_at_t0 temp1];
end
 
figure;
plot((0:(1/1000):0.999)*1000,abs(tv_tf_comp_at_t0));
title('Time Varying Transfer Function computed at the time instant t0=1us U19EC046')
xlabel('Frequency in Khz');
ylabel('Amplitude');
