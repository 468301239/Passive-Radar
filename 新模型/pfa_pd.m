clc;clear;close all;
k=1.380649*10^-23;
T=273.15+25;
N0=k*T;
a=N0;
SNR=[ 0.01 0.1 1];
p1=0:0.01:1;p0=1-p1;
N=8;
n=1:N;in=1;Tao=0;
wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n))); 
for i=1:length(SNR)
    d=SNR(i);
    for x=1:N
        for kx=x:N
            w=wn(kx);
            in=(in*besseli(0,2*a/N0*sqrt((N-x)./N*a^2+N0))./exp(d));
        end
        Tao=Tao+in;
        pd=(p1.*Tao)./(N.*p0+p1.*Tao);
    end
    plot(p1,pd,'LineWidth',1);
    hold on;
end
