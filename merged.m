%% image conversion
close all
img = imread('cameraman.tif');
% img = randi([1 100], 1000);
img_binary = dec2bin(img);
img_binary_reshaped = reshape(img_binary,[],1);
bits_per_symbol = 1;
extra_bits = bits_per_symbol - rem(numel(img_binary), bits_per_symbol)
extra = zeros(extra_bits, 1)
img_binary_reshaped = [img_binary_reshaped; extra];

payload = bin2dec(reshape(img_binary_reshaped, [numel(img_binary_reshaped)/bits_per_symbol, bits_per_symbol]));
M = 2^bits_per_symbol;

modulated = qammod(payload, M);
figure;
scatterplot(modulated)
figure;
noisy = awgn(modulated, 1);
% scatterplot(noisy)
demodulated = qamdemod(noisy, M);

img_binary_r = dec2bin(demodulated);
img_binary_r = reshape(img_binary_r, [numel(img_binary_r), 1]);
img_binary_r = img_binary_r(1:(numel(img_binary_r)-extra_bits), :);
img_binary_reshaped = reshape(img_binary_r, [numel(img_binary_r)/8, 8]);
img_dec_received = bin2dec(img_binary_reshaped);
img_dec_received_reshaped = reshape(img_dec_received, size(img));

imshow(uint8(img_dec_received_reshaped));



%% Multipath, transfer function and impulse response
clc;
clear all;
close all;
 
f = 1;
t = 0.01:1/100:1;
%Transmitted Signal.
x = cos(2*pi*f*t);
 
%J=No of Paths
J=2;
%B(j,t) is attenuation in jth path at t time.
B = rand(J,length(t));
%T(j,t) is time delay attenuation in jth path at t time.
T = rand(J,length(t));
%y is received multipath signal.
y = zeros(1,length(t));
 
%calculating received multipath signal.
for k=1:length(t)
    for j=1:J
        T(j,k)=T(j,k)*t(k);
        y(k)=y(k)+B(j,k)*exp(1i*2*pi*f*(t(k)-T(j,k)));
    end
end
 
%Plotting Transmitted Signal.
figure(1);
subplot(4,1,1);
plot(x);
title('U19EC002 Transmitted Signal');

%Plotting Received Signal.
subplot(4,1,2);
plot(real(y));
title('U19EC002 Received MultiPath Signal');
 
%Plotting Spectrum of Transmitted Signal.
subplot(4,1,3);
plot(abs(fft(x)));
title('U19EC002 Spectrum of Transmitted Signal');
 
%Plotting Spectrum of Received Signal.
subplot(4,1,4);
plot(abs(fft(real(y))));
title('U19EC002 Spectrum of Received MultiPath Signal');
 
%Calculating Transfer Function of Channel.
fs=100;
f=0:fs/101:(50*fs)/101;
Transfer_Function = zeros(length(t),length(f));
for k=1:length(t)
    for j=1:J
        Transfer_Function(k,:)=Transfer_Function(k,:)+B(j,k)*exp(-1i*2*pi*f*T(j,k));
    end
end
 
%Plotting Transfer Function for 42 time instants.
figure(2);
for k=1:42
    subplot(7,6,k);
    plot(real(Transfer_Function(k,:)));
    title(t(k));
end
 
%Plotting Impulse Response for 42 time instants.
figure(3);
for k=1:42
    subplot(7,6,k);
    %Impulse Response is inverse fourier transform of Transfer function.
    plot(real(ifft(Transfer_Function(k,:))));
    title(t(k));
end


%% QAM, QPSK, constellation diagram and BER vs SNR
clc
clear all;
close all;

% read the image
img=imread('cameraman.tif');

maxM = 6; % maximum symbol rate
maxSNR = 40; % maximum SNR value

% initializing an matrix to store BER for different symbol rates
BERQAM = zeros(maxM,maxSNR/2); 
BERPSK = zeros(maxM,maxSNR/2); 

QAMstr = 'QAM (M = %d)';
PSKstr = 'PSK (M = %d)';

% SNR ranges
t = 1:2:maxSNR;


figure(1);
subplot(maxM/3,3,1);
imshow(img);
title('Original');

figure(2);
subplot(maxM/3,3,1);
imshow(img);
title('Original');

% for loop to calculate stuff on different symbol rates
for m  = 2:maxM
   
% modulation order
ModOrd = 2^m;
symbolSize = log2(ModOrd);

% we need to pad zeros, why?
% suppose image size is 256*256 and each elment represents rgb value from 
% 0 to 255, i.e. 8 bits. Consider symbol rate 3, we need to represent image
% in array of 3 bits. if 256*256 is not divisible by 3 then we cannot
% reshape it into array of 3 bits. hence we need to pad the image matriz 
% so that number of elements in image size is divisible by symbol rate.

AddZero = rem(length(img),symbolSize);

if AddZero ~= 0
    img = [img; zeros(symbolSize - AddZero, numel(img)/length (img))];
end

% stuff meant to reshape the image matrix to array of m bits
binaryImage = de2bi(img);
reshapedImage = reshape(binaryImage, 8*length(binaryImage)/symbolSize, symbolSize);
img_dec = bi2de(reshapedImage);

% finnaly modulate
yQAM = qammod(img_dec, ModOrd);
yPSK = pskmod(double(img_dec), ModOrd);

% for each SNR value populate BER values
for s = 1:2:maxSNR
    nQAM = awgn(yQAM,s);
    nPSK = awgn(yPSK,s);
    zQAM = qamdemod(nQAM, ModOrd);
    zPSK = pskdemod(nPSK, ModOrd);
    [a,b] = biterr(img_dec,zQAM);
    % (s+1)/2 because SNR jumps 2 value each time 
    % and BER matrix has size [m, maxSNR]
    BERQAM(m,(s+1)/2) = 100*b;
    [c,d] = biterr(img_dec,zPSK);
    BERPSK(m,(s+1)/2) = 100*d;
end

% get the image matrix from demod, exatly opposite during modulation
QAM_dec = de2bi(zQAM);
QAM_rsp = reshape(QAM_dec, size(binaryImage));
QAMM = bi2de(QAM_rsp);
QAMM = uint8(reshape(QAMM,size(img)));

PSK_dec = uint8(de2bi(zPSK));
PSK_rsp = reshape(PSK_dec, size(binaryImage));
PSKK = bi2de(PSK_rsp);
PSKK = uint8(reshape(PSKK,size(img)));

% plot the scatter plot
scatterplot(nQAM);
scatterplot(nPSK);

figure(1);
subplot(maxM/3,3,m);
imshow(QAMM);
title(sprintf(QAMstr,ModOrd));

figure(2);
subplot(maxM/3,3,m);
imshow(PSKK);
title(sprintf(PSKstr,ModOrd));

figure(maxM+1);

subplot(212);
plot(t,BERPSK(m,:),'linewidth',2,'DisplayName',sprintf(PSKstr,ModOrd));
legend('m=2', 'm=3', 'm=4', 'm=5', 'm=6')
title('BER vs SNR (PSK) U19EC046');
xlabel('Signal to Noise Ratio');
ylabel('Bit Error Rate(%)');
hold on;

% subplot(211);
% plot(t,BERQAM((8-m),:),'linewidth',2,'DisplayName',sprintf(QAMstr,ModOrd));
% legend('m=2', 'm=3', 'm=4', 'm=5', 'm=6')
% title('BER vs SNR (QAM) U19EC046');
% xlabel('Signal to Noise Ratio');
% ylabel('Bit Error Rate(%)');
% hold on;


end

for m=2:maxM
    subplot(211);
    plot(t,BERQAM((8-m),:),'linewidth',2,'DisplayName',sprintf(QAMstr,ModOrd));
    axis([1 10 0 6])
    legend('m=2', 'm=3', 'm=4', 'm=5', 'm=6')
    title('BER vs SNR (QAM) U19EC046');
    xlabel('Signal to Noise Ratio');
    ylabel('Bit Error Rate(%)');
    hold on;
end


%% doppler spread in fast fading
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

%% delay spread in fast fading
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


%% doppler spread slow fading

clc;
clear all;
close all;
Tau0 =0;
f=1;
t = 0:0.01:100;
tx = cos(2*pi*f*t);
Tauj = [0.0042 0.0098 0.0030 0.0070];
BETA = [0.2691 0.4228 0.5479 0.9427];
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
xlim([1 1000]);
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
xlim([0 1000]);
title('U19EC002 Magnitude of Time varying transfer function');
subplot(2,1,2)
plot(phase(tf));
xlim([0 1000]);
title('U19EC002 Phase of Time varying transfer function');

%% delay spread slow fading


clear all;
close all;
Tau0=0;
t0=1;
Tauj = [0.0042 0.0098 0.0030 0.0070];
BETA = [0.2691 0.4228 0.5479 0.9427];
tv_tf_comp_at_t0=[];
for f=0:(1/1000):0.999
 temp=0;
 for i=1:length(Tauj)
 temp=temp+BETA(i)*exp(-1i*2*pi*f*Tau0)*exp(-1i*2*pi*f*Tauj(i)*t0);
 end
 tv_tf_comp_at_t0=[tv_tf_comp_at_t0 temp];
end
figure;
plot((0:(1/1000):0.999)*100000,abs(tv_tf_comp_at_t0));
title('U19EC002 Time Varying Transfer Function computed at the time instant t0=1us')
xlabel('Frequency in kHz');


%% DSSS and CDMA
clc;
clear all;
close all;
msg = randi([0,1],4,1);
msg_re = repmat(msg,[1,8]);
pn = generate(commsrc.pn('GenPoly', [3 2 0],'InitialStates',[0 0 1],'CurrentStates', [0 0 1],'Mask', [0 0 1], 'NumBitsOut',8));
rn = pn;
pn = repmat(pn.',[4,1]);
x_r = xor(msg_re,pn);
spread_i = x_r(:);
spread_i = uint8(spread_i);
snr = 4;
y1 = qammod(spread_i,4);
y2 = awgn(y1,snr);
y3 = qamdemod(y2,4);
y3 = reshape(y3,[],8);
despread_d = xor(y3,pn);
msg_rx = round(mean(despread_d,2));
BER = mean(abs(msg_rx - msg));
subplot(7,1,1);
plot(msg);
title('Message Signal (U19EC002)');
subplot(7,1,2);
plot(rn);
title('Random Noise (U19EC002)');
subplot(7,1,3);
plot(spread_i);
title('Spreaded Data (U19EC002)');
subplot(7,1,4);
plot(y1);
title('Modulated Signal (U19EC002)');
subplot(7,1,5);
plot(bi2de(y3));
title('Demodulated SIgnal (U19EC002)');
subplot(7,1,6);
plot(despread_d);
title('Despreaded Signal (U19EC002)');
subplot(7,1,7);
plot(msg_rx);
title('Received Signal (U19EC002)');

%BER vs SNR plot

msgbits = 1002;
colors = ['r', 'g', 'b'];
for j=1:3
msg = randi([0,1],msgbits,1);
msg_re = repmat(msg,[1,8]);
%PN sequence
if j==1
 H = commsrc.pn('Genpoly',[3 2 0],'InitialStates',[0 0 1],'CurrentStates',[0 0 1],'Mask',[0 0 1],'NumBitsOut',8);
 pn = generate(H);
elseif j==2
 H = commsrc.pn('Genpoly',[4 3 0],'InitialStates',[0 0 0 1],'CurrentStates',[0 0 0 1],'Mask',[0 0 0 1],'NumBitsOut',8);
 pn = generate(H);
else
 H = commsrc.pn('Genpoly',[5 3 0],'InitialStates',[0 0 0 0 1],'CurrentStates',[0 0 0 0 1],'Mask',[0 0 0 0 1],'NumBitsOut',8);
 pn = generate(H);
end
pn = reshape(pn,[1,8]);
pn = repmat(pn,msgbits,1);
x_r = xor(msg_re,pn);
spread_i = x_r(:);
spread_i = uint8(spread_i);
BER=[];
for snr = -10:1:10
 y1 = qammod(spread_i,4);
 y1 = awgn(y1,snr);
 y2 = qamdemod(y1,4);
 y2 = reshape(y2,[msgbits,8]);
 despread_d = xor(y2,pn);
 msg_rx = round(mean(despread_d,2));
 ber = mean(abs(msg_rx - msg));
 BER = [BER ber];
end
snr = -10:1:10;
plot(snr,BER,'DisplayName',sprintf('msg=%d',j),'color',colors(j));
title('BER v/s SNR plot (U19EC002)');
xlabel('SNR (dB)');
ylabel('BER');
grid on;
legend;
hold on;
end

%% DSSS part 2

% Code for BER vs SNR for 3 different message signals
clc;
clear all;
close all;
msgbits = 1002;
colors = ['r', 'g', 'b'];
for j=1:3
msg = randi([0,1],msgbits,1);
msg_re = repmat(msg,[1,8]);
%PN sequence
if j==1
 H = commsrc.pn('Genpoly',[3 2 0],'InitialStates',[0 0 1],'CurrentStates',[0 0 1],'Mask',[0 0 1],'NumBitsOut',8);
 pn = generate(H);
elseif j==2
 H = commsrc.pn('Genpoly',[4 3 0],'InitialStates',[0 0 0 1],'CurrentStates',[0 0 0 1],'Mask',[0 0 0 1],'NumBitsOut',8);
 pn = generate(H);
else
 H = commsrc.pn('Genpoly',[5 3 0],'InitialStates',[0 0 0 0 1],'CurrentStates',[0 0 0 0 1],'Mask',[0 0 0 0 1],'NumBitsOut',8);
 pn = generate(H);
end
pn = reshape(pn,[1,8]);
pn = repmat(pn,msgbits,1);
x_r = xor(msg_re,pn);
spread_i = x_r(:);
spread_i = uint8(spread_i);
BER=[];
for snr = -10:1:10
 y1 = qammod(spread_i,4);
 y1 = awgn(y1,snr);
 y2 = qamdemod(y1,4);
 y2 = reshape(y2,[msgbits,8]);
 despread_d = xor(y2,pn);
 msg_rx = round(mean(despread_d,2));
 ber = mean(abs(msg_rx - msg));
 BER = [BER ber];
end
snr = -10:1:10;
plot(snr,BER,'DisplayName',sprintf('msg=%d',j),'color',colors(j));
title('BER v/s SNR plot (U19EC002)');
xlabel('SNR (dB)');
ylabel('BER');
grid on;
legend;
hold on;
end

%% equalizer for bit stream
clc;
clear all;
close all;
M = 4; 
msg = randi([0 1], 5000, 1);
[n, m] = size(msg);
Mod=2;
data_vector = reshape(msg, [numel(msg)/Mod Mod]);
msg = bi2de(data_vector);
hMod = comm.QPSKModulator('PhaseOffset',0);
modmsg = step(hMod,msg); % Modulate using QPSK.
SNR= 40;
modmsg = awgn(modmsg , SNR);
scatterplot(modmsg);
legend('Modulated signal U19EC002');
trainlen = 200; % Length of training sequence
Tauj = [0.986 0.845 0.237 0.123];
Beta = [-0.1 0 -0.03 0.31];
chan = rayleighchan(1,0,Tauj, Beta);
chanCoeff = chan.AvgPathGaindB + 1i*chan.PathDelays;
filtmsg = filter(chanCoeff,1, modmsg); % Introduce channel distortion.
% Equalize the received signal.
eq1 = lineareq(8, lms(0.01)); % Create an equalizer object.
eq1.SigConst = step(hMod,(0:M-1)')'; % Set signal constellation.
[symbolest,yd] = equalize(eq1,filtmsg,modmsg(1:trainlen)); % Equalize.
scatterplot(filtmsg)
legend('Filtered signal U19EC002')
scatterplot(symbolest)
legend('Equilized signal U19EC002')
% Compute error rates with and without equalization.
hDemod = comm.QPSKDemodulator('PhaseOffset',0);
demodmsg_noeq = step(hDemod,filtmsg); % Demodulate unequalized signal.
demodmsg = step(hDemod,yd); % Demodulate detected signal from equalizer.
de = de2bi(demodmsg_noeq);
ber_without_eq = sum(mean(abs((de) - data_vector)))
de = de2bi(demodmsg);
ber_with_eq = sum(mean(abs((de) - data_vector)))
disp('Bit Error Rates with Equalisation')
disp(ber_with_eq)
disp('Bit Error Rates without Equalisation')
disp(ber_without_eq)

%% equalizer for image

clc;
clear all;
close all;
M = 4; 
id = imread('cameraman.tif');
figure('name','transmitted U19EC002');
imshow(id);
ida = id(:);
ib = de2bi(ida);
ib = ib(:);
msg=ib;
[n, m] = size(msg);
Mod=2;
data_vector = reshape(msg, [numel(msg)/Mod Mod]);
msg = bi2de(data_vector);
hMod = comm.QPSKModulator('PhaseOffset',0);
modmsg = step(hMod,msg); % Modulate using QPSK.
SNR= 40;
modmsg = awgn(modmsg , SNR);
scatterplot(modmsg);
legend('Modulated signal U19EC002');
trainlen = 200; % Length of training sequence
Tauj = [0.986 0.845 0.237 0.123];
Beta = [-0.1 0 -0.03 0.31];
chan = rayleighchan(1,0,Tauj, Beta);
chanCoeff = chan.AvgPathGaindB + 1i*chan.PathDelays;
filtmsg = filter(chanCoeff,1, modmsg); % Introduce channel distortion.
release(hMod)
% Equalize the received signal.
eq1 = lineareq(8, lms(0.01)); % Create an equalizer object.
eq1.SigConst = step(hMod,(0:M-1)')'; % Set signal constellation.
[symbolest,yd] = equalize(eq1,filtmsg,modmsg(1:trainlen)); % Equalize.
scatterplot(filtmsg)
legend('Filtered signal U19EC002')
scatterplot(symbolest)
legend('Equilized signal U19EC002')
% Compute error rates with and without equalization.
hDemod = comm.QPSKDemodulator('PhaseOffset',0);
demodmsg_noeq = step(hDemod,filtmsg); % Demodulate unequalized signal.
demodmsg = step(hDemod,yd); % Demodulate detected signal from equalizer.
% Demodulated image with using equalizer 
de = de2bi(demodmsg);
%Double to uint8 converion
r = uint8(de);
% r(length(r),:) = [];
r = r(:);
%Reshaping demodulated output
temp = r(1:length(r),:);
t = reshape(temp,[],8);
%Binary to Decimal conversion
k = bi2de(t);
%Reshaping into matrix of the size of image
p = reshape(k,256,256);
%Displaying output
figure('name','received image with equalizer U19EC002');
imshow(p);
ber_with_eq = sum(mean(abs(uint8(de) - data_vector)));
% Demodulated image without using equalizer 
de = de2bi(demodmsg_noeq);
%Double to uint8 converion
r = uint8(de);
% r(length(r),:) = [];
r = r(:);
%Reshaping demodulated output
temp = r(1:length(r),:);
t = reshape(temp,[],8);
%Binary to Decimal conversion
k = bi2de(t);
%Reshaping into matrix of the size of image
p = reshape(k,256,256);
%Displaying output
figure('name','received image without equalizer U19EC002');
imshow(p);
ber_without_eq = sum(mean(abs(uint8(de) - data_vector)));
disp('Bit error rates with and without equalizer:')
disp([ber_with_eq ber_without_eq])


