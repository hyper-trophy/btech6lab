

%===========================================================================%
% t = 0.01:0.01:1;
% f = 1;
% tx = cos(2*pi*f*t);
% J = 2
% B = rand(J,length(t));
% T = rand(J,length(t));
% rx = zeros(1, length(t));

% for i=1:length(t)
%     for j=1:J
%         T(j, i) = T(j, i)*t(i);
%         rx(i) = rx(i) + B(j, i)*exp(1i*2*pi*(t(i)-T(j, i)));
%     end
% end

% figure('name', 'test');
% subplot(4, 1, 1)
% plot(tx)

% subplot(4, 1, 2)
% plot(abs(rx))

% subplot(4, 1, 3)
% plot(abs(fft(tx)))

% subplot(4, 1, 4)
% plot(abs(fft(rx)))

% f = 0:1:100;
% tf = zeros(length(t), length(f));

% for i=1:length(t)
%     for j=1:J
%         tf(i, :) = tf(i, :) + B(j, i)*exp(-1i*2*pi*f*T(j, i));
%     end
% end

% plot(abs(ifft(tf(1,:))))

%=======================================================================================%
close all
img = imread('cameraman.tif');
img = randi([1 100], 1000);
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
