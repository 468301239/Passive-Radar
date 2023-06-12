clc;clear;close all;
k=1.380649*10^-23;T=273.15+25;
N=128;x0=N/2;n=0:N;N0=100*k*T;a=N0;SNR=3.16*a/N0;pi=3.14;
% wn=sqrt(N0/2).*(randn(1)+1i*randn(1));
% wn=sqrt(N0/2).*randn(N,1);


% % syms phi;syms x;
% fun=@(phi,x) exp(2*a/N0*abs(exp(-1i*phi)*wn))^(N-x+1);
% f=integral(fun,0,2*pi)/(2*pi);
% T=@(x) exp(-(N-x+1)*SNR)*f;
% Tao=integral(T,0,N)/N;

 T=@(x) exp(-(N-x+1).*SNR).*besseli(0,2*a/N0.*abs(wn)).^(N-x+1);
 Tao=integral(T,0,N)/N;


p1=0:0.01:1;
pfa=(p1*Tao)./(1-p1+p1*Tao);
plot(p1,pfa,'LineWidth',1);

% function[wn]=Wn(x)
% k=1.380649*10^-23;T=273.15+25;
% N0=100*k*T;
% wn=sqrt(N0/2)*(rand(fix(x))+1i*rand(fix(x)));
% end



% syms x
% ff=int( x*cos(x), 0, 1 );
% f=double(ff);
% y=1+f-1;
% plot(f,y)


