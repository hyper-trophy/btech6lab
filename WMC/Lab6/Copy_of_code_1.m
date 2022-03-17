% Set up parameters and signals.
M = 4; % Alphabet size for modulation

data = randi([0 1], 5000, 1);

[n, m] = size(data);

Mod=2;

data_vector = reshape(data, [numel(data)/Mod Mod]);

msg = bi2de(data_vector);

hMod = comm.QPSKModulator('PhaseOffset',0);
modmsg = step(hMod,msg); % Modulate using QPSK.
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
legend('Filtered signal')
scatterplot(symbolest)
legend('Equilized signal')

% Compute error rates with and without equalization.
hDemod = comm.QPSKDemodulator('PhaseOffset',0);
demodmsg_noeq = step(hDemod,filtmsg); % Demodulate unequalized signal.
demodmsg = step(hDemod,yd); % Demodulate detected signal from equalizer.
hErrorCalc = comm.ErrorRate; % ErrorRate calculator
ser_noEq = step(hErrorCalc, ...
    msg(trainlen+1:end), demodmsg_noeq(trainlen+1:end));
reset(hErrorCalc)
ser_Eq = step(hErrorCalc, msg(trainlen+1:end),demodmsg(trainlen+1:end));
disp('Symbol error rates with and without equalizer:')
disp([ser_Eq(1) ser_noEq(1)])