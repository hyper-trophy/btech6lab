%% doppler spread
clc;
clear all;
close all;
%beta is assumed as the one that is not varying with time
% nop is the no of paths
%txsignal is the transmit signal
% rx signal is the received signal
%tvtf is the time varying transfer function computed at the f0 
%frequency
%TAUJ is the rate at which the delay is changing
%The f0 is in MHz and sampling time is in micro seconds
Tau0 =0;
f=1;
t = 0:0.01:100;
tx = cos(2*pi*f*t);
Tauj = [0.62 1.84 0.86 0.37];
BETA = [0.23 0.17 0.23 0.44];
f_shift=abs(-f*Tauj+f);
rx =[];
tf = [];
for t =0:0.01:100
 temp =0;
 temp1 =0;
 for i=1:length(Tauj)
 temp = temp+ BETA(i)*cos(2*pi*f_shift(i)*t);
 temp1 = temp1+BETA(i)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(i)*t);
 end
 rx = [rx temp];
 tf = [tf temp1];
end
figure(1)
subplot(2,2,1)
plot(tx)
axis([1 1000 -1.5 1.5]);
title('U19EC002 Transmitted signal');
xlabel('Samples with Ts=1/100 micro seconds');
subplot(2,2,2)
plot(rx)
axis([1 1000 -1.5 1.5]);
title('U19EC002 Received signal after subject to multipath');
xlabel('Samples with Ts=1/100 micro seconds');
fre = (0:1:length(rx)-1)/100;
subplot(2,2,3)
plot(fre, abs(fft(tx)));
xlim([0 2]);
title('U19EC002 Spectrum of Transfer signal');
xlabel('Frequency in MHz');
subplot(2,2,4)
plot(fre, abs(fft(rx)));
xlim([0 2]);
title('U19EC002 Spectrum of Received signal');
xlabel('Frequency in MHz');
figure(2)
subplot(2,1,1)
plot(abs(tf));
axis([0 1000 0 1.5]);
title('U19EC002 Magnitude of Time varying transfer function');
subplot(2,1,2)
plot(phase(tf));
axis([0 1000 -25 0]);
title('U19EC002 Phase of Time varying transfer function');

%% delay spread

clc;
clear all;
close all;
Tau0=0;
t0=1;
Tauj = [0.62 1.84 0.86 0.37];
BETA = [0.23 0.17 0.23 0.44];
tv_tf_comp_at_t0=[];
for f=0:(1/1000):0.999
 temp=0;
 for i=1:length(Tauj)
 temp=temp+BETA(i)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(i)*t0);
 end
 tv_tf_comp_at_t0=[tv_tf_comp_at_t0 temp];
end
figure;
plot((0:(1/1000):0.999)*1000,abs(tv_tf_comp_at_t0));
title('U19EC002 Time Varying Transfer Function computed at the time instant t0=1us')
xlabel('Frequency in kHz');
