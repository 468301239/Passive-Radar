%%I(Z;X)
clc;clear;close all;
PI=3.14;
deltax=0.05; 
for N=[128,1024]
x=0:deltax:N-deltax;
SNR =-30:70;%信噪比（dB）
alpha=1;
rho2 = 10.^(SNR/10);%信噪比(数字形式)
I=zeros(1,length(SNR));
for i=1:length(SNR)
    sum=0;
    sum_p_w_x=0;
    snr=rho2(i);
    N0=2*alpha^2/snr;
%     sum_p_w_x=-1/(2*alpha^2)*(-erf(snr*sqrt(N/(4*snr+2)))...
%                               +exp(N/2)*erf(sqrt(N/2))...
%                               -exp(N/2)*erf((1+snr)*sqrt(N/(4*snr+2))));
    for kx=1:length(x)
        sum_p_w_x=sum_p_w_x+...
                  1/sqrt(2*PI*(N*N0^2+2*(N-x(kx))*alpha^2*N0))*deltax;
    end

    for kx=1:length(x)
         p_x_w=1/sum_p_w_x/sqrt(2*PI*(N*N0^2+2*(N-x(kx))*alpha^2*N0));
         sum=sum-p_x_w*log2(p_x_w+eps)*deltax;
    end
    I(1,i)=log2(N)-sum;
end
plot(SNR,I,'LineWidth',2);
hold on
end
legend('N=128','N=1024','N=30000');
xlabel('SNR');ylabel('I(Z;X)');
title('不同时间带宽积恒模散射目标的距离信息与SNR的关系');