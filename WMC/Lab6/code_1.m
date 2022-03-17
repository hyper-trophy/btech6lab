clc
clear all

data = randi([0 1], 1000, 1);

[n m] = size(data);

M=2;

data_vector = reshape(data, [numel(data)/M M]);

decimal_data = bi2de(data_vector);

transmitted = pskmod(decimal_data, 2^M);

% hMod = comm.QPSKModulator('PhaseOffset',0);
% transmitted = step(hMod,decimal_data); % Modulate using QPSK.

trainlen = 500; % Length of training sequence

scatterplot(transmitted);

SNR = 10;

noisy = awgn(transmitted, SNR);

scatterplot(noisy);

Tauj = [0.62 1.84 0.86 0.37]*1e-6;

BETA = [1 1 1 1];

chan = rayleighchan(1,0,[0.1 0.2 0.2 0.1],[-0.5 0 -1 -0.4])
% chan = [.986; .845; .237; .123+.31i]/10; % Channel coefficients

received = filter(chan, noisy);

scatterplot(received);

% alg = lms(0.01)

% eqobj = lineareq(8,alg);

% y = equalize(eqobj, received, transmitted(1:200));
eq1 = lineareq(4, lms(0.01)); % Create an equalizer object.
% eq1.SigConst = step(hMod,(0:M-1)')'; % Set signal constellation.
[symbolest,yd] = equalize(eq1,received,transmitted); % Equalize.


scatterplot(symbolest);

% % Set up parameters and signals.
% M = 4; % Alphabet size for modulation
% msg = randi([0 M-1],1500,1); % Random message
% hMod = comm.QPSKModulator('PhaseOffset',0);
% modmsg = step(hMod,msg); % Modulate using QPSK.
% trainlen = 500; % Length of training sequence
% chan = [.986; .845; .237; .123+.31i]; % Channel coefficients
% filtmsg = filter(chan,1,modmsg); % Introduce channel distortion.

% % Equalize the received signal.
% eq1 = lineareq(8, lms(0.01)); % Create an equalizer object.
% eq1.SigConst = step(hMod,(0:M-1)')'; % Set signal constellation.
% [symbolest,yd] = equalize(eq1,filtmsg,modmsg(1:trainlen)); % Equalize.

% % Plot signals.
% h = scatterplot(filtmsg,1,trainlen,'bx'); hold on;
% scatterplot(symbolest,1,trainlen,'g.',h);
% scatterplot(eq1.SigConst,1,0,'k*',h);
% legend('Filtered signal','Equalized signal',...
%    'Ideal signal constellation');
% hold off;

% % Compute error rates with and without equalization.
% hDemod = comm.QPSKDemodulator('PhaseOffset',0);
% demodmsg_noeq = step(hDemod,filtmsg); % Demodulate unequalized signal.
% demodmsg = step(hDemod,yd); % Demodulate detected signal from equalizer.
% hErrorCalc = comm.ErrorRate; % ErrorRate calculator
% ser_noEq = step(hErrorCalc, ...
%     msg(trainlen+1:end), demodmsg_noeq(trainlen+1:end));
% reset(hErrorCalc)
% ser_Eq = step(hErrorCalc, msg(trainlen+1:end),demodmsg(trainlen+1:end));
% disp('Symbol error rates with and without equalizer:')
% disp([ser_Eq(1) ser_noEq(1)])