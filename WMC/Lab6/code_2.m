clear all;
close all;
% Set up parameters and signals.
img=imread('cameraman.tif');

M = 4; % Alphabet size for modulation
symbolSize=2;

AddZero = rem(length(img),symbolSize);

if AddZero ~= 0
    img = [img; zeros(symbolSize - AddZero, numel(img)/length (img))];
end

% stuff meant to reshape the image matrix to array of m bits
binaryImage = de2bi(img);
reshapedImage = reshape(binaryImage, 8*length(binaryImage)/symbolSize, symbolSize);
msg = bi2de(reshapedImage);


% hMod = comm.QPSKModulator('PhaseOffset',0);
% modmsg = step(hMod,msg); % Modulate using QPSK.

modmsg = pskmod(double(msg), M);

trainlen = length(msg)/2; % Length of training sequence

Tauj = [0.986 0.845 0.237 0.123];
Beta = [-0.1 0 -0.03 0.31];

chan = rayleighchan(1,0,Tauj, Beta);
chanCoeff = chan.AvgPathGaindB + 1i*chan.PathDelays;
filtmsg = filter(chanCoeff,1, modmsg); % Introduce channel distortion.

% Equalize the received signal.
eq1 = lineareq(8, lms(0.01)); % Create an equalizer object.
% eq1.SigConst = step(hMod,(0:M-1)')'; % Set signal constellation.
[symbolest,yd] = equalize(eq1,filtmsg,modmsg(1:trainlen)); % Equalize.


% scatterplot(filtmsg)
% legend('Filtered signal')
% scatterplot(symbolest)
% legend('Equilized signal')

% Compute error rates with and without equalization.
% hDemod = comm.QPSKDemodulator('PhaseOffset',0);
% demodmsg_noeq = step(hDemod,filtmsg); % Demodulate unequalized signal.
% demodmsg = step(hDemod,yd); % Demodulate detected signal from equalizer.

demodmsg = pskdemod(yd, M);

QAM_dec = de2bi(demodmsg);
QAM_rsp = reshape(QAM_dec, size(binaryImage));
QAMM = bi2de(QAM_rsp);
QAMM = uint8(reshape(QAMM,size(img)));

% QAM_dec2 = de2bi(demodmsg_noeq);
% QAM_rsp2 = reshape(QAM_dec2, size(binaryImage));
% QAMM2 = bi2de(QAM_rsp);
% QAMM2 = uint8(reshape(QAMM,size(img)));

figure;
imshow(img)
figure;
imshow(QAMM)
% hErrorCalc = comm.ErrorRate; % ErrorRate calculator
% ser_noEq = step(hErrorCalc, ...
%     msg(trainlen+1:end), demodmsg_noeq(trainlen+1:end));
% reset(hErrorCalc)
% ser_Eq = step(hErrorCalc, msg(trainlen+1:end),demodmsg(trainlen+1:end));
% disp('Symbol error rates with and without equalizer:')
% disp([ser_Eq(1) ser_noEq(1)])