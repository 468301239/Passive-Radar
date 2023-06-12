clear;close all;clc;
N = 16;N0 = 1;
SNR =5;

P_FA = 0:0.01:1;
%K = length(P_FA);
T0 = erfinv(1-2*P_FA);
P_D = 0.5*erfc((T0/sqrt(N0*(1+2*SNR))-N*SNR/sqrt(2*N*(1+2*SNR))));
plot(P_FA,P_D,'LineWidth',1);hold on;

xlabel('Ðé¾¯¸ÅÂÊ');
ylabel('¼ì²â¸ÅÂÊ');
hold on


