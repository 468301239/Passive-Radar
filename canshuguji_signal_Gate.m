%% %%%%%%距离信息%%%%%%%%%%%%%%%%%%
clc;clear;close all;
N=16;
deltax = 0.01;
x=-N/2:deltax:N/2-deltax;
x0=0;
n=-N/2:1:N/2-1;
SNR = -10:30;
rho2 = 10.^(SNR/10);
num =1000;
alpha=1;
I1=zeros(1,length(SNR));
I2=zeros(1,length(SNR));

parfor i=1:length(SNR)
        N0=2*alpha^2/rho2(i);
        H1=zeros(1,num);
        H2=zeros(1,num);
        for run=1:num
            p_w_x1=zeros(1,length(x));
            p_w_x2=zeros(1,length(x));
            wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
            for ii=1:length(x)
                p_w_x1(ii)=besseli(0,2*alpha/N0*abs(alpha*gate(x(ii)-x0)+sum(wn.*gate(n-x(ii)))));
                p_w_x2(ii)=besseli(0,2*alpha/N0*abs(alpha*sinc(x(ii)-x0)+sum(wn.*sinc(n-x(ii)))));
            end
            p_x_w1=p_w_x1./sum(p_w_x1*deltax);
            p_x_w2=p_w_x2./sum(p_w_x2*deltax);
            H1(run)=sum(-p_x_w1.*log2(p_x_w1)*deltax);
            H2(run)=sum(-p_x_w2.*log2(p_x_w2)*deltax);
        end
        I1(i)=log2(N)-mean(H1);
        I2(i)=log2(N)-mean(H2);
        %     Imax(i)=log2(N)-0.5*log2(2*pi*exp(1)*CRB);
end
plot(SNR,smooth(I1),'LineWidth',1);hold on;grid on;
plot(SNR,smooth(I2),'LineWidth',1);
% title('门宽=0.05');
% plot(SNR,Imax,'LineWidth',1);

%% %%%%%%%%%先验后验%%%%%%%%%%%%
clc;clear;close all;
N=16;
deltax = 0.01;
x=-N/2:deltax:N/2-1;
x0=0;
n=-N/2:1:N/2-1;
SNR = 5;
rho2 = 10.^(SNR/10);
num =1;
alpha=1;
I1=zeros(1,length(SNR));
I2=zeros(1,length(SNR));

for i=1:length(SNR)
    N0=2*alpha^2/rho2(i);
    H1=zeros(1,num);
    H2=zeros(1,num);
    for run=1:num
        p_w_x1=zeros(1,length(x));
        p_w_x2=zeros(1,length(x));
        wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
        for ii=1:length(x)
            p_w_x1(ii)=besseli(0,2*alpha/N0*abs(alpha*gate(x(ii)-x0)+sum(wn.*gate(n-x(ii)))));
            p_w_x2(ii)=besseli(0,2*alpha/N0*abs(alpha*sinc(x(ii)-x0)+sum(wn.*sinc(n-x(ii)))));
        end
        p_x_w1=p_w_x1./sum(p_w_x1*deltax);
        p_x_w2=p_w_x2./sum(p_w_x2*deltax);
    end
end
plot(x,p_x_w1*deltax,'LineWidth',1);hold on;
plot(x,p_x_w2*deltax,'LineWidth',1);hold on;
legend('gate函数','sinc函数');
title('5dB下后验PDF');


%% gate*gate
clc;clear;close all;

N=16;
deltax =1;
x=-N/2:deltax:N/2-1;
x0=0;
n=-N/2:N/2-1;
y1=gate(n-x);
y2=gate(n-x0);
y3=(gate(n-x).*gate(n-x0));
y4=gate(x-x0);
figure
subplot(211)
% plot(x,y1);hold on
% plot(x,y2);hold on
plot(x,y3);hold on;legend('gate(n-x_0)*gate(n-x)')
subplot(212)
plot(x,y4);hold on;legend('gate(x-x_0)')
% legend('gate(n-x)','gate(n-x0)','gate(n-x)*gate(n-x0)');


%% gate
clc;clear;close all;

N=8;
deltax =0.01;
x=-N/2:deltax:N/2;
y1=gate(x);
y2=sinc(x);
plot(x,y1,x,y2);
legend('gate','sinc')