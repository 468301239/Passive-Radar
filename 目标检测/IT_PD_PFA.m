%%%IT下的PD和PFA
clear ;close all;clc;
N=64;N0=1;P_FA=0:0.00001:1;
x=0:0.5:2;e=2.718;
snr=[0.01 0.316 1];%-10db -5db 0db
for i=1:length(snr)
    d=snr(i);
    up=P_FA.*exp(((1+d).*N.*d^2)./(1+2.*d));  
    down=1+P_FA.*(exp(((1+d).*N.*d^2)./(1+2.*d))-1);
    P_D=up./down;
    plot(P_FA,P_D,'LineWidth',1);
    set(gca,'XScale','log');%对数坐标Pfa
    xlabel('虚警概率');
    ylabel('检测概率');
    legend('SNR=-10dB','SNR=-5dB','SNR=0dB');
    hold on;
end  

% pfa=0:0.00001:1;  %虚警概率 
% % N=100;
% d=0.316;
%  M=[32,64,128,256];% 时间带宽积N
% %M=[0.01,0.1,0.3,1,0];%dB 低信噪比
% for ii=1:length(M)
%     N=M(ii);
%      c=exp((N.*d^2.*(1+d))./(1+2*d));
%      c1=1+(c-1).*pfa;
%      pd=(c.*pfa)./c1;
%      hold on
%     plot(pfa,pd);
%     set(gca,'XScale','log');%对数坐标Pfa
%     title('N 取不同值时的ROC曲线');
% end
% grid on; xlabel('PFA'); ylabel('PD');
% legend('d=-20dB','d=-10dB','d=-5dB','d=0dB','SNR=0')%低
% legend('N=32','N=64','N=128','N=256')%
