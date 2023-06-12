clc;clear;close all;
k=1.380649*10^-23;T=273.15+25;
N=1024;x0=N/2;n=0:N;N0=100*k*T;a=1;SNR=1;
wn=sqrt(N0/2)*(randn(1)+1i*randn(1));
% wn=sqrt(1/SNR)*(randn(1)+1i*randn(1));
% wn=0.57;

% %%%x<x0
% x=0:0.01:x0;
% y1=x.*(SNR^2-1-N0*SNR^4);
% plot(x,y1,'r','LineWidth',2);
% hold on;
% %%%x>x0
% x=x0:0.01:N;
% y2=x.*(SNR^2-1-(N0+1)*SNR^4);
% plot(x,y2,'b','LineWidth',2);



%x<x0
x=1:0.01:x0;
P1=SNR.*x+(x0-x+1).*log(besseli(0,2.*SNR.*abs(wn./a)))+...
          (N-x0+1).*log(besseli(0,2.*SNR.*abs(1+wn./a)));
% semilogy(x,P1,'r','LineWidth',2);
plot(x,P1,'r');
hold on


%x>x0
x2=x0:0.01:N;
P2=SNR.*x2+(N-x2+1).*log(besseli(0,2*SNR*abs(1+wn./a)));
% semilogy(x2,P2,'b','LineWidth',2);
plot(x2,P2,'b');
