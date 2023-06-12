clc;clear;close all;
SNR=10;N=16;N0=1;a=1;x0=3*N/4;x=0:0.001:x0;n=ceil(x):N;x2=x0:0.001:N;
wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));


in=2.*SNR.*(abs(N-x0+sum(wn)./a));
P_up=exp(SNR.*x).*besseli(0,in);
plot(x,P_up);
hold on;
in2=2.*SNR.*(abs(N-x2+sum(wn)./a));
P_up=exp(SNR.*x2).*besseli(0,in2);
plot(x2,P_up);

% for j=1:length(SNRS)
%     figure(j);
%     SNR=SNRS(j);
%     in=2.*SNR.*(abs(abs(x-x0)+sum(wn)./a));
% %     in=2.*SNR.*(abs(abs(x-x0)+1));
%     P_up=exp(SNR.*x).*besseli(0,in);
% %     subplot(2,3,i);
%     plot(x,P_up);
%     hold on;
% end


%x0<x
% for j=1:length(SNRS)
%     figure(j);
%     SNR=SNRS(j);
%     for i=1:length(X0)
%         x0=X0(i);
%         if x0<=x
%             in=2.*SNR.*(abs(N-x+sum(wn)./a));
%         else
%             in=2.*SNR.*(abs(N-x0+sum(wn)./a));
%         end
%         P_up=exp(SNR.*x).*besseli(0,in);
%         subplot(2,3,i);
%         plot(x,P_up);
%         hold on;
%     end
% end


% if x0<=x
%     in=2.*SNR.*(abs(N-x+sum(wn)./a));
% elseif x0>=x
%     in=2.*SNR.*(abs(N-x0+sum(wn)./a));
% end
%%%x0<=x

