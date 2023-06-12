%%%IT和NP的ROC  PD-PFA
clear all;close all;clc;
N=128;P_FA=0:0.000001:1;
snr=[0.1 0.3162];
%%%NP
for i=1:length(snr)
    d=snr(i);
    P_D_NP=0.5*erfc((sqrt(2*N)*erfcinv(2*P_FA)-d*N)/(sqrt(2*N+4*N*d)));   
    plot(P_FA,P_D_NP,'--','LineWidth',1);
    set(gca,'XScale','log');%对数坐标Pfa
    hold on
end


%%%IT
for i=1:length(snr)
    d=snr(i);
    up=P_FA.*exp(((1+d).*N.*d^2)./(1+2.*d));  
    down=1+P_FA.*(exp(((1+d).*N.*d^2)./(1+2.*d))-1);
    P_D_IT=up./down;
    plot(P_FA,P_D_IT,'LineWidth',1);
    set(gca,'XScale','log');%对数坐标Pfa
    hold on
end  

xlabel('PFA');ylabel('PD');
legend('SNR=-10dB   NP','SNR=-5dB   NP','SNR=-10dB  IT','SNR=-5dB   IT');
