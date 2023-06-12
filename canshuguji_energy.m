%%I(Z;X) EE CRB
clc;clear;close all;
deltax=0.01; 
num=100;
SNR =-15:30;%信噪比（dB）
alpha=1;
rho2 = 10.^(SNR/10);%信噪比(数字形式)
h=zeros(1,num);
I=zeros(1,length(SNR));
I_max=zeros(1,length(SNR));
EE=zeros(1,length(SNR));
CRB=zeros(1,length(SNR));
MSE=zeros(1,length(SNR));
row=zeros(1,num);
summ=zeros(1,num);

for N=[128]
    x0=N/2;
    x=0:deltax:N-deltax;
    n = 0 : 1 : N-1;
    p_w_x = zeros(1,length(x)); 
    for i=1:length(SNR)
        snr=rho2(i);
        N0=2*alpha^2/snr;
        for run=1:num
            for kx=1:length(x)
                p_w_x(1,kx)=1/sqrt((N+2*(N-x(kx))*snr))*...
                          exp(-(x(kx)-x0)^2*snr^2/2/(N+2*(N-x(kx))*snr));
            end
            p_x_w = p_w_x./sum(p_w_x*deltax); 
            h(1,run)=-sum(p_x_w.*log2(p_x_w+eps))*deltax; 

            maxp=max(p_w_x);
            row(1,run)=find(p_w_x==maxp);
            Location=(row-1)*deltax;        
            summ(1,run)=(Location(1,run)-x0)^2; 
        end
        CRB(1,i)=N/snr;%-(N+2*(N-x0)*snr)^2/((2-N)*snr^2+2*snr^3*(N-x0));
        I(1,i)=log2(N)-mean(h);
        I_max(1,i)=log2(N)-0.5*log2(2*pi*exp(1)*CRB(1,i));
        EE(1,i)=2^(2*mean(h))/(2*pi*exp(1));   
        MSE(1,i)=sum(summ)/num;
    end
%     plot(SNR,I,'LineWidth',2);hold on;grid on
%     plot(SNR,I_max,'--','LineWidth',2);hold on

    semilogy(SNR,EE,'LineWidth',2);hold on;grid on
    semilogy(SNR,CRB,'--','LineWidth',2);hold on
    semilogy(SNR,MSE,'--','LineWidth',2);hold on
end
legend('EE','CRB','MSE');
xlabel('SNR');ylabel('方差');
title('恒模散射目标归一化时延方差');

% legend('N=128距离信息','N=128距离信息上界','N=1024距离信息','N=1024距离信息上界');
% xlabel('SNR');ylabel('I(Z;X)');
% title('不同时间带宽积恒模散射目标的距离信息与SNR的关系');



%%EE
% clc;clear;close all;
% PI=3.14;e=2.718;
% deltax=0.01; 
% for N=[128]
% x=0:deltax:N;
% SNR =-30:50;%信噪比（dB）
% alpha=1;
% rho2 = 10.^(SNR/10);%信噪比(数字形式)
% I=zeros(1,length(SNR));
% for i=1:length(SNR)
%     h_x_z=0;
%     sum_p_w_x=0;
%     snr=rho2(i);
%     N0=2*alpha^2/snr;
% %     sum_p_w_x=-1/(2*alpha^2)*(-erf(snr*sqrt(N/(4*snr+2)))...
% %                               +exp(N/2)*erf(sqrt(N/2))...
% %                               -exp(N/2)*erf((1+snr)*sqrt(N/(4*snr+2))));
%     for kx=1:length(x)
%         sum_p_w_x=sum_p_w_x+...
%                   1/sqrt(2*PI*(N*N0^2+2*(N-x(kx))*alpha^2*N0))*...
%                   exp(-(N-x(kx))^2*alpha^4/2/(N*N0^2+2*(N-x(kx))*alpha^2*N0))*deltax;
%     end
% 
%     for kx=1:length(x)
%          p_x_w=1/sum_p_w_x/sqrt(2*PI*(N*N0^2+2*(N-x(kx))*alpha^2*N0))*...
%                 exp(-(N-x(kx))^2*alpha^4/2/(N*N0^2+2*(N-x(kx))*alpha^2*N0));
%          h_x_z=h_x_z-p_x_w*log2(p_x_w+eps)*deltax;
%     end
%     SIGMAEE(1,i)=2^(2*h_x_z)/(2*PI*e);
% 
% %     for kx=1:length(x)
% %         f_x=1/sqrt(2*PI*(N*N0^2+2*(N-x(kx))*alpha^2*N0))*...
% %                   exp(-(N-x(kx))^2*alpha^4/2/(N*N0^2+2*(N-x(kx))*alpha^2*N0));
% %         f_x_plusdeltax=1/sqrt(2*PI*(N*N0^2+2*(N-x(kx)-deltax)*alpha^2*N0))*...
% %                   exp(-(N-x(kx)-deltax)^2*alpha^4/2/(N*N0^2+2*(N-x(kx)-deltax)*alpha^2*N0));
% %         f_x_minusdeltax=1/sqrt(2*PI*(N*N0^2+2*(N-x(kx)+deltax)*alpha^2*N0))*...
% %                   exp(-(N-x(kx)+deltax)^2*alpha^4/2/(N*N0^2+2*(N-x(kx)+deltax)*alpha^2*N0));
% %         f_2dao(1,kx)=f_x_plusdeltax+f_x_minusdeltax-2*f_x;
% %     end
% %     FI_x=-mean(f_2dao);
% 
% %     down=2*N+4*snr*(N-x);
% %     FI_x=-mean(2*snr^2./(N+2*snr*(N-x)).^2-2*snr^2./down...
% %        +16*snr^3*(N-x)./down.^2-32*snr^4*(N-x).^2./down.^3);
% 
%     down=N+2*snr*(N-x);
%     FI_x=-mean(2*snr^2./down.^2-snr^2./down+4*snr^3*(N-x)./down.^2-...
%                 8*snr^4*(N-x).^2./down.^3);
% 
%     CRB(1,i)=1/FI_x;
% end
% semilogy(SNR,SIGMAEE,'LineWidth',2);
% hold on
% semilogy(SNR,CRB,'LineWidth',2);
% hold on
% end
% legend('EE','CRB');
% xlabel('SNR');ylabel('方差');


% %%CRB
% clc;clear;close all;
% PI=3.14;e=2.718;
% deltax=0.05; 
% for N=[128,1024]
% x=0:deltax:N-deltax;
% SNR =-30:50;%信噪比（dB）
% alpha=1;
% rho2 = 10.^(SNR/10);%信噪比(数字形式)
% I=zeros(1,length(SNR));
% for i=1:length(SNR)
%     snr=rho2(i);
%     N0=alpha^2/snr;
% %     syms x
% %     f(x)=1/sqrt(2*PI*(N*N0^2+2*(N-x)*alpha^2*N0))*...
% %             exp(-(N-x)^2*alpha^4/2/(N*N0^2+2*(N-x)*alpha^2*N0));
% %     CRB(1,i)=-1/diff(f(x),2);
% 
%     for kx=1:length(x)
%         f_x=1/sqrt(2*PI*(N*N0^2+2*(N-x(kx))*alpha^2*N0))*...
%                   exp(-(N-x(kx))^2*alpha^4/2/(N*N0^2+2*(N-x(kx))*alpha^2*N0));
%         f_x_plusdeltax=1/sqrt(2*PI*(N*N0^2+2*(N-x(kx)-deltax)*alpha^2*N0))*...
%                   exp(-(N-x(kx)-deltax)^2*alpha^4/2/(N*N0^2+2*(N-x(kx)-deltax)*alpha^2*N0));
%         f_x_minusdeltax=1/sqrt(2*PI*(N*N0^2+2*(N-x(kx)+deltax)*alpha^2*N0))*...
%                   exp(-(N-x(kx)+deltax)^2*alpha^4/2/(N*N0^2+2*(N-x(kx)+deltax)*alpha^2*N0));
%         f_2dao(1,kx)=f_x_plusdeltax+f_x_minusdeltax-2*f_x;
%     end
%     FI_x=-mean(f_2dao);
% 
% % down=2*N+4*snr*(N-x);
% % FI_x=-mean(2*snr^2./(N+2*snr*(N-x)).^2-2*snr^2./down...
% %        +16*snr^3*(N-x)./down.^2-32*snr^4*(N-x).^2./down.^3);
% 
% CRB(1,i)=1/FI_x;
% end
% semilogy(SNR,CRB,'LineWidth',2);
% hold on
% end
