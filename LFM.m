%%%%%%%%% 线性调频信号LFM的时频特性和模糊特性 %%%%%%%%

% @Author:37yuany

% @Time:2021/4/21

%% 1.LFM时域波形

clc; clear; close all;

tic

f0 = 0; %雷达中心频率

T = 1e-7;%脉宽

B = 1e9; %带宽

fs = 3*B;

Ts = 1/fs;

N = T/Ts;

k = B/T;

t = linspace(-T/2,T/2,N);

y = 1*exp(1j*(2*pi*f0*t + pi*k*t.^2));

figure;

plot(t*1e6,real(y));xlabel('时间（us）');ylabel('幅度');

title('LFM信号时域波形（实部）');

grid on; axis tight;

%% 2.LFM频谱图

Sf = fftshift(fft(y));

f = linspace(-fs/2,fs/2,N);

figure(2);

plot(f*1e-6,abs(Sf)./max(max(abs(Sf))));

xlabel('频率（MHz）')

ylabel('归一化幅度频谱');

title('LFM信号的频谱图');

grid on;axis tight;

%% 3.LFM的模糊函数

Grid = 1000;

t = -T:T/Grid:T;

f = -B:B/Grid:B;

[tau,fd]=meshgrid(t,f);

var1=T-abs(tau);

var2=pi*(fd-k*tau).*var1;

var2=var2+eps;

amf=abs(sin(var2)./var2.*var1/T);

amf=amf/max(max(amf));

var3=pi*k*tau.*var1;

taul=abs(sin(var3)./var3.*var1);

taul=taul/max(max(taul));

mul=T.*abs(sin(pi*fd.*T)./(pi*fd.*T));

mul=mul/max(max(mul));

figure

mesh(tau.*1e6,fd*1e-6,amf);

xlabel ('时延（us）');ylabel ('多普勒频率（MHz）');zlabel ('归一化幅度');

title('LFM信号三维模糊函数');

grid on;axis tight;

figure

contour(tau.*1e6,fd*1e-6,amf);

xlabel ('时延（us）');ylabel ('多普勒频率（MHz）');

title('LFM信号等高图');

grid on;axis tight;

figure

plot(fd.*1e-6,mul(:,Grid+1))

xlabel ('多普勒频率（MHz）');ylabel ('归一化幅度');

title('LFM信号速度切面图');

grid;axis tight;

figure

plot(t.*1e6,taul(Grid+1,:));

xlabel ('时延（us）');ylabel ('归一化幅度');

title('LFM信号距离切面图');

grid on;axis tight;

toc
