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
