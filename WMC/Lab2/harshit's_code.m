clc;
clear all;
close all;
 
f = 1;
t = 0.01:1/100:1;
%Transmitted Signal.
x = cos(2*pi*f*t);
 
%J=No of Paths
J=2;
%B(j,t) is attenuation in jth path at t time.
B = rand(J,length(t));
%T(j,t) is time delay attenuation in jth path at t time.
T = rand(J,length(t));
%y is received multipath signal.
y = zeros(1,length(t));
 
%calculating received multipath signal.
for k=1:length(t)
    for j=1:J
        T(j,k)=T(j,k)*t(k);
        y(k)=y(k)+B(j,k)*exp(1i*2*pi*f*(t(k)-T(j,k)));
    end
end
 
%Plotting Transmitted Signal.
figure(1);
subplot(4,1,1);
plot(x);
title('U19EC002 Transmitted Signal');

%Plotting Received Signal.
subplot(4,1,2);
plot(real(y));
title('U19EC002 Received MultiPath Signal');
 
%Plotting Spectrum of Transmitted Signal.
subplot(4,1,3);
plot(abs(fft(x)));
title('U19EC002 Spectrum of Transmitted Signal');
 
%Plotting Spectrum of Received Signal.
subplot(4,1,4);
plot(abs(fft(real(y))));
title('U19EC002 Spectrum of Received MultiPath Signal');
 
%Calculating Transfer Function of Channel.
fs=100;
f=0:fs/101:(50*fs)/101;
Transfer_Function = zeros(length(t),length(f));
for k=1:length(t)
    for j=1:J
        Transfer_Function(k,:)=Transfer_Function(k,:)+B(j,k)*exp(-1i*2*pi*f*T(j,k));
    end
end
 
%Plotting Transfer Function for 42 time instants.
figure(2);
for k=1:42
    subplot(7,6,k);
    plot(real(Transfer_Function(k,:)));
    title(t(k));
end
 
%Plotting Impulse Response for 42 time instants.
figure(3);
for k=1:42
    subplot(7,6,k);
    %Impulse Response is inverse fourier transform of Transfer function.
    plot(real(ifft(Transfer_Function(k,:))));
    title(t(k));
end

%second code
% clc;
% clear all;
% close all;
% f=1;
% number_of_path=2;
% received_signal=[];
% t=0:1/100:1;
% transmitted_signal=cos(2*pi*f*t);
% z=1;
% for t=0:1/100:1
%  temp=0;
%  for p=1:1:number_of_path
%  beta(p)=rand;
%  delay(p)=rand*t;
%  temp=temp+beta(p)*exp(1i*2*pi*f*(t-delay(p)));
%  end
%  BETACOL{z}=beta;
%  DELAYCOL{z}=delay;
%  beta=0;
%  delay=0;
%  received_signal=[received_signal temp];
%  z=z+1;
% end;
% figure(1);
% subplot(4,1,1);
% plot(transmitted_signal);
% title('transmitted signal U19EC023');
% xlabel('time');
% ylabel('amplitude');
% subplot(4,1,2);
% plot(real(received_signal));
% title('received signal after multipath U19EC023');
% xlabel('time');
% ylabel('amplitude');
% subplot(4,1,3);
% plot(abs(fft(transmitted_signal)));
% title('spectrum of transmitted signal U19EC023');
% xlabel('frequency');
% ylabel('amplitude');
% subplot(4,1,4);
% plot(abs(fft(received_signal)));
% title('spectrum of received signal after multipath U19EC023');
% xlabel('frequency');
% ylabel('amplitude');
% fs=100;
% u=1;
% for f=0:fs/101:(50*fs)/101
%  received_signal=[];
%  temp=0;
%  z=1;
%  for t=0:1/100:1
%  temp=0;
%  for p=1:1:number_of_path
%  temp=temp+BETACOL{z}(p)*exp(1i*2*pi*f*(t-DELAYCOL{z}(p)));
%  end
%  received_signal=[received_signal temp];
%  z=z+1;
%  end
 
%  t=0:1/100:1;
%  timevaryingTF_at_freq_f{u}=received_signal.*exp(-1i*2*pi*f*t);
%  u=u+1;
% end
% TEMP=cell2mat(timevaryingTF_at_freq_f');
% for i=1:1:101
%  u=TEMP(:,i);
%  u1=[u;transpose(u(length(u):-1:2)')];
%  timevaringIR_at_time_t{i} = ifft(u1);
% end
% Transfer_Function_MATRIX=abs(cell2mat(timevaryingTF_at_freq_f'));
% Impulse_Response_MATRIX=cell2mat(timevaringIR_at_time_t);
% s=[2:2:8];
% for i=1:1:4
%  figure(2)
%  subplot (2,2,i)
%  plot(Impulse_Response_MATRIX(1:1:101,s(i)))
%  title(strcat('(U19EC023) t=',num2str((s(i)-1)/100)))
%  ylabel('amplitude');
%  xlabel('frequency');
 
%  figure(3)
%  subplot(2,2,i)
%  plot(Transfer_Function_MATRIX(:,s(i)))
%  title(strcat('(U19EC023) t=',num2str((s(i)-1)/100)))
%  ylabel('amplitude');
%  xlabel('frequency');
% end