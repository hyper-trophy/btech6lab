close all
img = imread('cameraman.tif');
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
img_dec_received_reshaped = reshape(img_dec_received, [256, 256]);

imshow(uint8(img_dec_received_reshaped));
