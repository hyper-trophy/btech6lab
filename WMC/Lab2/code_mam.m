clc
clear all

f=5;
nop=2; % initial value of nop = 2
received_signal=[];
t=0:1/100:1; % choosing sampling freg=100Hz
transmitted_signal=cos(2*pi*f*t); % transmitted signal
z=1;
for t=0:1/100:1
    temp=0;
    for p=1:1:nop
        beta(p)=rand/5; % for every delayed singal there will be 10
        delay(p)=rand*t/5; % delay of each multipth component generated
        temp=temp+beta(p)*exp(1i*2*pi*f*(t-delay(p)));
    end
BETACOL{z}=beta;
DELAYCOL{z}=delay;
beta=0;
delay=0;
received_signal=[received_signal temp];
z=z+1;
end

save CONSTANTS BETACOL DELAYCOL

figure('name', 'U19EC046')

subplot(4,1,1)

plot(transmitted_signal)

title('U19EC046 Transmitted Signal');

subplot(4,1,2)

plot(real(received_signal))

title('U19EC046 received signal after multipath');
subplot(4,1,3)

plot(abs(fft(transmitted_signal)))

title('U19EC046 spectrum of the transmitted signal');
subplot(4,1,4)

plot(abs(fft(real(received_signal))))

title('U19EC046 spectrum of the received signal after multipath');

hold

load CONSTANTS
fs=100;
u=1;
for f=0:fs/101:(50*fs)/101
    received_signal=[];
    temp=0; 
    z=1;
    for t=0:1/100:1
        temp=0;
        for p=1:1:nop
            temp=temp+BETACOL{z}(p)*exp(1i*2*pi*f*(t-DELAYCOL{z}(p)));
        end
    received_signal=[received_signal temp];
    z=z+1;
    end
    % The impulse response of the time-varying channel is computed as folld
    t=0:1/100:1;
    timevaryingTF_at_freq_f{u}=received_signal.*exp(-1i*2*pi*f*t);
    u=u+1;
end

TEMP=cell2mat(timevaryingTF_at_freq_f');
for i=1:1:101
    u=TEMP(:,i);
    u1=[u;transpose(u(length(u):-1:2)')];
    timevaringIR_at_time_t{i} = ifft(u1);
end 

TFMATRIX=abs(cell2mat(timevaryingTF_at_freq_f'));
IRMATRIX=cell2mat(timevaringIR_at_time_t);
s=[2:2:8];

figure('name', 'U19EC046 Impulse Response');
for i=1:1:4
    subplot (2,2,i)
    plot(IRMATRIX(1:1:101,s(i)))
    title(strcat('t=',num2str((s(i)-1)/100)))
end

figure('name', 'U19EC046 Transfer Function');
for j=1:1:4
    subplot (2,2,j)
    plot(real(ifft(IRMATRIX(1:1:101,s(j)))))
    title(strcat('t=',num2str((s(j)-1)/100)))
end

%% Raghav's code
clc;
clear all;
close all;

f = 1;
t = 0.01:1/100:1;
x = cos(2*pi*f*t);
J=3;
B = rand(J,length(t));
T = rand(J,length(t));
y = zeros(1,length(t));

for k=1:length(t)
    for j=1:J
        T(j,k)=T(j,k)*t(k);
        y(k)=y(k)+B(j,k)*exp(1i*2*pi*f*(t(k)-T(j,k)));
    end
end

figure(1);
subplot(4,1,1);
plot(x);
title('U19EC043 Transmitted Signal');

subplot(4,1,2);
plot(real(y));
title('U19EC043 Received MultiPath Signal');

subplot(4,1,3);
plot(abs(fft(x)));
title('U19EC043 Spectrum of Transmitted Signal');

subplot(4,1,4);
plot(abs(fft(real(y))));
title('U19EC043 Spectrum of Received MultiPath Signal');

fs=100;
f=0:fs/101:(50*fs)/101;
Transfer_Function = zeros(length(t),length(f));
for k=1:length(t)
    for j=1:J
        Transfer_Function(k,:)=Transfer_Function(k,:)+B(j,k)*exp(-1i*2*pi*f*T(j,k));
    end
end

figure(2);
% Change this part according to you(How many values of k?)
for k=1:4
    subplot(2,2,k);
    plot(real(Transfer_Function(k,:)));
    title(t(k));
end

figure(3);
% Change this part according to you(How many values of k?)
for k=1:4
    subplot(2,2,k);
    plot(real(ifft(Transfer_Function(k,:))));
    title(t(k));
end

    