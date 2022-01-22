clc;
clear all;
close all;
M = 64;
img = imread('cameraman.tif');
figure('name', 'transmitted')
imshow(img);
a = img(:);
b = de2bi(a);
c = double(b);
d = c(:);
y = pskmod(d,M);
e = pskdemod(y,M);
f = reshape(e,[65536,8]);
g = uint8(f);
h = bi2de(g);
x = reshape(h,256,256);
figure('name', 'received')
imshow(x);
scatterplot(y);
scatterplot(e);