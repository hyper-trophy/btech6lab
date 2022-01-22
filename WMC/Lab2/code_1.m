% mera
clc;
clear all;

fc = 1/100;
dt = 1/100;
t = 0:dt:100;
carrier = cos(2*pi*fc*t);

beta = rand(1, 10);

fshift = rand(1, 10)/100;

received = zeros(1, length(t));

for j=1:length(beta)
    received = received + beta(j)*cos(2*pi*(fc + (-1)^randi(2)*fshift(j))*t);
end

carrierFft = fft(carrier);

receivedFft = fft(received);

subplot(4, 1, 1)
plot(t, carrier)

subplot(4, 1, 2)
plot(t, carrierFft)
axis([-5 120 -500 10000])

subplot(4, 1, 3)
plot(t, received)

subplot(4, 1, 4)
plot(t, receivedFft)
axis([-5 120 -500 10000])
