clc;clear;close all;
pi=3.1416;
alpha=1;
N=2^4;
n=-N/2:N/2-1;
SNR=-20:1:25;
rho=10.^(SNR/10);
deltax=0.01;
x=-N/2:deltax:N/2-0.01;
Num=100;
x0=0;
I_snr=zeros(1,length(SNR));
I_snrr=zeros(1,length(SNR));
sigma_CRB=zeros(1,length(SNR));
for snr=1:length(SNR)
    N0=2*alpha^2/rho(snr);
    I_run=zeros(1,Num);
    for run=1:Num
        p1=zeros(1,length(x));
        wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
        for kx=1:length(x)
             p1(kx)=2*pi*besseli(0,2*alpha/N0*abs(alpha*sinc(x(kx))+sum(wn.*sinc(n-x(kx)))));
        end
        p2=p1./sum(p1*deltax);
        I_run(run)=sum(-p2.*log2(p2)*deltax);
    end
    sigma_CRB(snr) = 3/(pi^2*rho(snr));
    I_snr(snr)=log2(N)-mean(I_run);%2^(2*mean(I_run))/(2*pi*e);
    I_snrr(snr) = log2(N) - 1/2*log2(2*pi*exp(1)*sigma_CRB(snr));
end
hold on;plot(SNR,I_snr,'-','LineWidth',1);grid on;
%  hold on;plot(SNR,I_snrr,'k--','LineWidth',1);grid on;
xlabel('SNR(dB)');ylabel('I(Z;X)(bit)');
legend('距离互信息','距离互信息上界')