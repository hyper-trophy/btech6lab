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