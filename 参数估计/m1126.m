clc;clear;close all;
SNR=1;N=16;N0=1;a=1;x0=N/2;
n=0:N;
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
wn=0.57;

%x<x0
x=0:0.01:x0;
P1=SNR.*x+(x0-x).*log(besseli(0,2.*SNR.*abs(wn/a)))+...
          (N-x0).*log(besseli(0,2.*SNR.*abs(1+wn/a)));
semilogy(x,P1,'r');
% plot(x,P1,'r');
hold on

%x>x0
x2=x0:0.01:N;
P2=SNR.*x2+(N-x2).*log(besseli(0,2*SNR*abs(1+wn/a)));
semilogy(x2,P2,'b');
% plot(x2,P2,'b');












% %x<x0
% x=0:0.01:x0;
% in1=2.*SNR.*(abs(sum(wn)/a));
% in2=2.*SNR.*(abs(1+sum(wn)/a));
% P_up1=exp(SNR.*x).*((besseli(0,in1)).^(x0-x)).*((besseli(0,in2)).^(N-x0));
% semilogy(x,P_up1,'r');
% % plot(x,P_up1);
% hold on
% 
% randn(20)
% %x>x0
% x=x0:0.01:N;
% n=ceil(x):N;
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% in=2.*SNR.*(abs(1+sum(wn)/a));
% P_up2=exp(SNR.*x).*((besseli(0,in)).^(N-x));
% semilogy(x,P_up2,'b');






% for x=0:0.01:x0
%     n1_1=ceil(x):x0;
%     wn1_1=sqrt(N0/2)*(randn(1,length(n1_1))+1i*randn(1,length(n1_1)));
%     in1_1=2.*SNR.*(abs(sum(wn1_1)./a));
%     I1_1=I1_1.*(besseli(0,in1_1));
%     
%     n1_2=ceil(x0):N;
%     wn1_2=sqrt(N0/2)*(randn(1,length(n1_2))+1i*randn(1,length(n1_2)));
%     in1_2=2.*SNR.*(abs(1+sum(wn1_2)./a));
%     I1_2=I1_2.*(besseli(0,in1_2));
%     pup_1=exp(SNR.*x).*I1_1.*I1_2;
%     plot(x,pup_1);
%     hold on;
% end
