clc;clear;close all;
k=1.380649*10^-23;T=273.15+25;
N=128;x0=N/2;n=0:N;N0=100*k*T;a=N0;SNR=a/N0;
wn=sqrt(N0/2)*(randn(1)+1i*randn(1));

%x<x0
x=0:0.01:x0;
P1=(exp(x.*SNR^2).*(1+SNR^2).^(-x).*SNR^2)./(exp(x0.*SNR^2)-1);
% semilogy(x,P1,'r','LineWidth',2);
% subplot(1,2,1);
plot(x,P1,'r','LineWidth',2);
hold on

%x>x0
x=x0:0.01:N;
% P2=(exp(x.*SNR^2).*(1+SNR^4*(abs(1+wn/a))^2).^(-x)*(SNR^2-log(1+SNR^4*(abs(1+wn/a))^2)))./...
%     (exp(x0.*SNR^2)*(1+SNR^4*(abs(1+wn/a))^2).^(-x0)-1);
P2=(exp(x.*SNR^2).*(1+SNR^4*(abs(1+wn/a))^2).^(-x)*(SNR^2))./...
    (exp(x0.*SNR^2)-1);


% semilogy(x,P2,'b','LineWidth',2);
% subplot(1,2,2);
plot(x,P2,'b--','LineWidth',2);