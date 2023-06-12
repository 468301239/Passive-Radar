clear all;close all;clc;
N=4;P_FA=0:0.0001:1;
% snr=[0 1 1.59 3.16 10];
% SNR=1:100;
% snr=10*log10(SNR);
snr=1:10;
for i=1:length(snr)
    P_D=0.5*erfc((sqrt(2*N)*erfinv(1-2*P_FA)-snr(i)*N)/(sqrt(2*N+4*N*snr(i))));  
    plot(P_FA,P_D,'LineWidth',1);%set(gca,'XScale','log');%对数坐标Pfa
    hold on;
end
xlabel('虚警概率');
ylabel('检测概率');
% legend('SNR=0','SNR=3','SNR=10','SNR=13');