clear;
close all;
clc;

N =4;
N0 = 1;
%P_FA = 0:0.0001:1;
P_FA = 0:0.000001:1;


snr = [1 10^0.5];
for i=1:2
    SNR = snr(i);
    T0 = erfinv(1-2*P_FA);
    X = N*SNR/sqrt(2*N*(1+2*SNR));
    P_D = 0.5*erfc(T0/sqrt(1+2.*SNR)-X);
    %NP检测的ROC
    figure(1)
    semilogx(P_FA,P_D,'--')
    xlabel('虚警概率');
    ylabel('检测概率');
    axis([1e-6 1 0 1]);
    hold on
    legend('SNR=0dB','SNR=5dB')

end

%信息论的ROC
for i= 1:2
    SNR = snr(i);
    P1 = exp(0.5*N*SNR^2./(1+2*SNR))./(exp(0.5*N*SNR^2./(1+2*SNR))*sqrt(1+2*SNR)+1./P_FA+1);
    P_D_I = P1*exp(0.5*N*SNR^2)/((1-P1)*sqrt(1+2*SNR)+P1*exp(0.5*N*SNR^2));
    semilogx(P_FA,P_D_I,'*')
    hold on;
end
