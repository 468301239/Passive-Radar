%% 后验的分子绘图
clc;clear;close all;
k=1.380649*10^-23;T=273.15+25;
N=30000;x0=N/2;n=1:N;N0=100*k*T;SNR=0.1;
a=sqrt(SNR*N0);
wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
w=mean(wn);
% x<x0分子形状
x=0:0.1:x0;
P1=SNR.*x+...
    (x0-x+1).*log(besseli(0,2*SNR*abs(w/a)))+...
    (N-x0+1).*log(besseli(0,2*SNR*abs(1+w/a)));
% semilogy(x,P1,'r','LineWidth',2);
% subplot(1,2,1);
plot(x,P1,'r','LineWidth',2);
hold on
grid on

% x>x0分子形状
x=x0:0.1:N;
P2=SNR.*x+(N-x+1).*log(besseli(0,2*SNR*(1+w/a)));
% semilogy(x,P2,'b','LineWidth',2);
% subplot(1,2,2);
plot(x,P2,'b','LineWidth',2);
% ylim([1.3, 1.9]);

%% 低SNR 后验近似 约分母
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;
% N=16;x0=N/2;n=1:2048;N0=100*k*T;SNR=1;
% a=sqrt(SNR*N0);
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% w=mean(wn);
% 
% % x<x0
% in1=@(x) exp(SNR.*x).*(1+SNR).^(-x);
% down1=integral(in1,0,x0);
% x=0:0.1:x0;
% p_left=(exp(SNR.*x).*(1+SNR).^(-x))./down1;
% plot(x,p_left,'r','LineWidth',2);
% hold on
% % x>x0
% in2=@(x) exp(SNR.*x).*(1+(SNR*abs(1+w/a))^2).^(-x);
% down2=integral(in2,x0,N);
% x=x0:0.1:N;
% p_right=exp(SNR.*x).*(1+(SNR*abs(1+w/a))^2).^(-x)./down2;
% plot(x,p_right,'r','LineWidth',2);

%% 低SNR 后验未近似
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;
% N=32;x0=N/2;n=1:N;N0=k*T;SNR=1;
% a=sqrt(SNR*N0);
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% w=sum(wn);
% 
% % x<x0
% in1=@(x) exp(SNR.*x).*(besseli(0,2*SNR*abs(w/a))).^(x0-x+1).*(besseli(0,2*SNR*abs(1+w/a))).^(N-x0+1);
% down1=integral(in1,0,x0);
% x=0:0.01:x0;
% p_left=exp(SNR.*x).*(besseli(0,2*SNR*abs(w/a))).^(x0-x+1).*(besseli(0,2*SNR*abs(1+w/a))).^(N-x0+1)./down1;
% plot(x,p_left,'r','LineWidth',2);
% hold on
% % x>x0
% in2=@(x) exp(SNR.*x).*(besseli(0,2*SNR*abs(1+w/a))).^(N-x+1);
% down2=integral(in2,x0,N);
% x=x0:0.01:N;
% p_right=exp(SNR.*x).*(besseli(0,2*SNR*abs(1+w/a))).^(N-x+1)./down2;
% plot(x,p_right,'b','LineWidth',2);
