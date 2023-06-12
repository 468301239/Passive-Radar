clc;clear;close all;
k=1.380649*10^-23;
T=273.15+25;
N0=k*T;
a=N0;p1=0:0.01:1;
SNR=[0.01 0.1 1 10];
pfa=0:0.01:1;
for i=1:length(SNR)
    d=SNR(i);
    pd=marcumq(sqrt(2*d),sqrt(-2*log(pfa)));
    plot(pfa,pd,'LineWidth',1);
    hold on
end
xlabel('P_F_A');ylabel('P_D');
legend( 'SNR=-20dB', 'SNR=-10dB','SNR=   0dB','SNR= 10dB');
set(gca,'FontName','Times New Roman','FontSize',12)