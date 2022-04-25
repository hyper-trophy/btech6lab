clc;
clear all;
close all;
M = input('Enter modulation order: ');
%Reading image
id = imread('cameraman.tif');
figure('name','transmitted ()');
imshow(id);
%Input matrix to 1d row vector
ida = id(:);
%Decimal to Binary conversion
ib1 = de2bi(ida);
ibd = double(ib1);
%Reshaping column matrix
ib = ibd(:);
%Modulation
x = mod(length(ib),log2(M));
no_of_zeros=log2(M)-x;
ib5 = [ib.'(zeros(1,no_of_zeros)];
%Reshaping column matrix
ib2 = reshape(ib5,[],log2(M));
ib3 = bi2de(ib2);
%Modulation
pm = pskmod(ib3,M);
%Demodulation
scatterplot(pm);

pdm = pskdemod(pm,M);

%Reshaping demodulated output
de = de2bi(pdm);
%Double to uint8 converion
r = uint8(de);
r = r(:);
%Reshaping demodulated output
temp = r(1:length(r)-no_of_zeros,:);
t = reshape(temp,[],8);

%Binary to Decimal conversion
k = bi2de(t);
%Reshaping into matrix of the size of image
p = reshape(k,256,256);
%Displaying output
figure('name','received ');
imshow(p);