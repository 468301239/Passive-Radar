clc;clear;close all;
SNRS=[0.1 1 10 100];N=8;N0=1;a=1;x=0:1:N;n=ceil(x):N;X0=0:2:8;

for j=1:length(SNRS)
    figure(j);
    SNR=SNRS(j);
    for i=1:length(X0)
        x0=X0(i);
        P_up=exp(-abs(x-x0).*(SNR+2*a));
        subplot(2,3,i);
        plot(x,P_up);
        hold on;
    end
end