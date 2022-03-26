clc;
clear all;
close all;
 
M = 4;
m = log2(M);
%Reading image
id = imread('cameraman.tif');
figure(1)
subplot(1,3,1)
imshow(id);
title('U19EC002 Input');
%Decimal to Binary conversion
ib = de2bi(id);
ibd = double(ib);
ib = ibd(:);
x = mod(length(ib),log2(M));
zo = log2(M)-x;
ib5 = [ib.' zeros(1,zo)];
%Reshaping column matrix
ib2 = reshape(ib5,[],log2(M));
ib3 = bi2de(ib2);
%Modulation
Mod = comm.QPSKModulator('PhaseOffset',0);
qm = pskmod(ib3,M);
 
snr = 20;
o_noise = awgn(qm,snr);
 
Tau = 0.48;
PdB = -2;
chan = rayleighchan(1,0,Tau,PdB);
chanCoeff = chan.AvgPathGaindB + 1i*chan.PathDelays;
filtmsg = filter(chanCoeff,1,o_noise);
trainlen = 500;
 
% Equalize the received signal.
eq1 = lineareq(8, lms(0.01)); % Create an equalizer object.
eq1.SigConst = Mod((0:M-1)')'; % Set signal constellation.
[symbolest,yd] = equalize(eq1,filtmsg,o_noise(1:trainlen)); % Equalize.
 
h = scatterplot(filtmsg,1,trainlen,'bx'); 
hold on;
scatterplot(symbolest,1,trainlen,'g.',h);
legend('Filtered signal','Equalized signal');
hold off;
 
%Without Equalizer
z_wo_eq = pskdemod(filtmsg,M);
de = de2bi(z_wo_eq,m);
[a,b] = biterr(ib2,de);
ber = b;
r = uint8(de);
r = r(:);
%Reshaping demodulated output
temp = r(1:length(r)-zo,:);
t = reshape(temp,[],8);
%Binary to Decimal conversion
k = bi2de(t);
%Reshaping into matrix of the size of image
out = reshape(k,256,256);
figure(1)
subplot(1,3,2)
imshow(out);
title('Output without equalization');
 
%With Equalizer
z_w_eq = pskdemod(yd,M);   
de1 = de2bi(z_w_eq,m);
[a1,b1] = biterr(ib2,de1);
ber1 = b1;
r = uint8(de1);
r = r(:);
%Reshaping demodulated output
temp = r(1:length(r)-zo,:);
t = reshape(temp,[],8);
%Binary to Decimal conversion
k = bi2de(t);
%Reshaping into matrix of the size of image
out = reshape(k,256,256);
figure(1)
subplot(1,3,3)
imshow(out);
title('Output with equalization');
 
sprintf('BER with and without equalizer are %f and %f respectively',ber1,ber)