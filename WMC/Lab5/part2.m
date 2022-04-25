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