clc
clear all;
close all;
img=imread('cameraman.tif');
maxM = 6;
maxSNR = 40;
BERQAM = zeros(maxM,maxSNR/2);
BERPSK = zeros(maxM,maxSNR/2);
QAMstr = 'QAM (M = %d)';
PSKstr = 'PSK (M = %d)';
t = 1:2:maxSNR;
figure(1);
subplot(maxM/3,3,1);
imshow(img);
title('Original U19EC002');
figure(2);
subplot(maxM/3,3,1);
imshow(img);
title('Original U19EC002');
color = ['b','g','r','c','m','y']
for m = 2:maxM
 ModOrd = 2^m;
 SymSize = log2(ModOrd);
 AddZero = rem(length(img),SymSize);
 if AddZero ~= 0;
 img = [img; zeros(SymSize - AddZero, prod(size(img))/length(img))];
 end
 img_bin = de2bi(img);
 img_rsp = reshape(img_bin, 8*length(img_bin)/SymSize, SymSize);
 img_dec = bi2de(img_rsp);
 yQAM = qammod(img_dec, ModOrd);
 yPSK = pskmod(double(img_dec), ModOrd);
 for s = 1:2:maxSNR
 nQAM = awgn(yQAM,s);
 nPSK = awgn(yPSK,s);
 zQAM = qamdemod(nQAM, ModOrd);
 zPSK = pskdemod(nPSK, ModOrd);
 [a,b] = biterr(img_dec,zQAM);
 BERQAM(m,(s+1)/2) = 100*b;
 [c,d] = biterr(img_dec,zPSK);
 BERPSK(m,(s+1)/2) = 100*d;
 end
 QAM_dec = de2bi(zQAM);
 QAM_rsp = reshape(QAM_dec, size(img_bin));
 QAMM = bi2de(QAM_rsp);
 QAMM = uint8(reshape(QAMM,size(img)));
 PSK_dec = uint8(de2bi(zPSK));
 PSK_rsp = reshape(PSK_dec, size(img_bin));
 PSKK = bi2de(PSK_rsp);
 PSKK = uint8(reshape(PSKK,size(img)));
 
 scatterplot(nQAM);
 title('U19EC002');
 scatterplot(nPSK);
 title('U19EC002');
 figure(1);
 subplot(maxM/3,3,m);
 imshow(QAMM);
 title(sprintf(QAMstr,ModOrd));
 figure(2);
 subplot(maxM/3,3,m);
 imshow(PSKK);
 title(sprintf(PSKstr,ModOrd));
 figure(maxM+1);
 subplot(211);
 plot(t,BERQAM(m,:),color(m-1),'linewidth',2,'DisplayName',sprintf(QAMstr,ModOrd));
 title('BER vs SNR (QAM) U19EC002');
 xlabel('Signal to Noise Ratio');
 ylabel('Bit Error Rate(%)');
 legend;
 hold on;
 subplot(212);
 plot(t,BERPSK(m,:),color(m-1),'linewidth',2,'DisplayName',sprintf(PSKstr,ModOrd));
 title('BER vs SNR (PSK) U19EC002');
 xlabel('Signal to Noise Ratio');
 ylabel('Bit Error Rate(%)');
 legend;
 hold on;
end