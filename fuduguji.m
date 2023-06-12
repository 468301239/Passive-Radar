%% CRB EE MSE
clc;clear;close all;
alpha=1;
num=500;N=16;deltax = 0.05;x0=10;
SNR=-10:15;
n=0:1:N-1;
rho2 = 10.^(SNR/10);%信噪比
x=0:deltax:N-1;
x1=0:deltax:x0-deltax;
x2=x0:deltax:N-1;

p_z_x=zeros(1,length(x));
h=zeros(1,num);
I=zeros(1,length(SNR));
EE=zeros(1,length(SNR));
for i=1:length(SNR)
     N0=2*alpha^2/rho2(i);
     summ1=0;
     for run=1:num
        wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
        for ii=1:length(x1)
%             wn1=sqrt(N0/2)*(randn(1,N-round(ii*deltax))+1i*randn(1,N-round(ii*deltax)));
%             p_z_x(1,ii)=exp(rho2(i)/2*x(ii)).*...
%                         besseli(0,rho2(i)*abs(N-x0+1/alpha*sum(wn1)));
            p_z_x(1,ii)=exp(rho2(i)/2*x(ii)).*...
                        besseli(0,rho2(i)*abs(N-x0+1/alpha*sum(wn(1,round(x(ii)+1:N)))));
        end
        for jj=(length(x1)+1):length(x)
%             wn2=sqrt(N0/2)*(randn(1,N-round(jj*deltax))+1i*randn(1,N-round(jj*deltax)));
%             p_z_x(1,jj)=exp(rho2(i)/2*x(jj)).*...
%                         besseli(0,rho2(i)*abs(N-x(jj)+1/alpha*sum(wn2)));
            p_z_x(1,jj)=exp(rho2(i)/2*x(jj)).*...
                        besseli(0,rho2(i)*abs(N-x(jj)+1/alpha*sum(wn(1,round(x(jj)+1:N)))));
        end
        p_x_z=p_z_x/sum(p_z_x*deltax);
        h(1,run)=-sum(p_x_z.*log2(p_x_z+eps))*deltax;
        
        maxp1=max(p_x_z);
        row1(1,run)=find(p_x_z==maxp1);
        Location1=(row1-1).*deltax;        
        summ1=summ1+(Location1(1,run)-x0)^2; 

     end
     I(1,i)=log2(N)-mean(h);
     EE(1,i)=2^(2*mean(h))/(2*pi*exp(1));
     MSE(1,i)=summ1/num;
     CRB(1,i)=1/((0.5*rho2(i)/2+1.5*rho2(i)/2).^2);
end
% semilogy(SNR, EE,'LineWidth',1);hold on;grid on
% semilogy(SNR, CRB,'--','LineWidth',1);hold on;grid on
% semilogy(SNR, MSE,'--','LineWidth',1);hold on;grid on
plot(SNR, I,'LineWidth',2);hold on;grid on







%% 距离信息xh
clc;clear;close all;
num=100;N=32;
x0=0;
SNR=-15:15;
alpha=1;
rho2 = 10.^(SNR/10);%信噪比
deltax = 0.1;
x=-N/2:deltax:N/2-1;
x1=-N/2:deltax:x0-deltax;
x2=x0:deltax:N/2-1;
step1=0;
step2=0;

p_z_x_left=zeros(1,length(x1));
p_z_x_right=zeros(1,length(x2));
h_left=zeros(1,num);
h_right=zeros(1,num);
I=zeros(1,length(SNR));
p_B_z_x_left=zeros(1,length(x1));
p_B_z_x_right=zeros(1,length(x2));
I_B=zeros(1,length(SNR));
EE=zeros(1,length(SNR));
CRB=zeros(1,length(SNR));

for i=1:length(SNR)
     N0=2*alpha^2/rho2(i);
     for run=1:num
        for ii=1:length(x1)
            wn1=sqrt(N0/2)*(randn(1,N/2-ii)+1i*randn(1,N/2-ii));
            p_z_x_left(1,ii)=exp(rho2(i)/2*x1(ii)-step1).*...
                             besseli(0,rho2(i)*abs(N-x0+1/alpha*sum(wn1))-step2);
%             p_B_z_x_left(1,ii)=exp(rho2(i)/2*x1(ii)-step1)*...
%                                besseli(0,rho2(i)*abs(N)-step2);
        end
        p_x_z_left=p_z_x_left/sum(p_z_x_left*deltax);
        for jj=1:length(x2)
            wn2=sqrt(N0/2)*(randn(1,N/2-jj)+1i*randn(1,N/2-jj));
            p_z_x_right(1,jj)=exp(rho2(i)/2*x2(jj)-step1).*...
                              besseli(0,rho2(i)*abs(N-x2(jj)+1/alpha*sum(wn2))-step2);
%             p_B_z_x_right(1,jj)=exp(rho2(i)/2*x2(jj)-step1)*...
%                               besseli(0,rho2(i)*abs(N)-step2);
        end
        p_x_z_right=p_z_x_right/sum(p_z_x_right*deltax);

%         plot(x1, p_x_z_left,'LineWidth',1);hold on
%         plot(x2, p_x_z_right,'LineWidth',1);hold on
%         p_B_x_z_left=p_B_z_x_left/sum(p_B_z_x_left*deltax);
%         p_B_x_z_right=p_B_z_x_right/sum(p_B_z_x_right*deltax);
%         h_B=-sum(p_B_x_z_left.*log2(p_B_x_z_left+eps)*deltax)...
%             -sum(p_B_x_z_right.*log2(p_B_x_z_right+eps)*deltax);

%         h_left(1,run)=-sum(p_x_z_left.*log2(p_x_z_left+eps))*deltax;
%         h_right(1,run)=-sum(p_x_z_right.*log2(p_x_z_right+eps))*deltax;
%         h=h_left+h_right;
p_x_z=p_z_x/sum(p_z_x*deltax);
          h(1,run)=-sum(p_x_z.*log2(p_x_z+eps))*deltax;
     end
     I(1,i)=log2(N)-mean(h);
     EE(1,i)=2^(2*mean(h))/(2*pi*exp(1));
     CRB(1,i)=1/(0.5*((rho2(i)/2)+(3*rho2(i)/2)));
%      I_B(1,i)=log2(N)-mean(h_B);
%     CRB(1,i)=
end
semilogy(SNR, EE,'LineWidth',1);hold on;grid on
semilogy(SNR, CRB,'LineWidth',1);hold on;grid on
% plot(SNR, I,'LineWidth',1);hold on;grid on
% plot(SNR, I_B,'--','LineWidth',1);hold on





% %% 后验
% clc;clear;close all;
% num=100;N=32;
% x0=N/2;n=0:N-1;
% SNR=20;
% alpha=1;
% rho2 = 10.^(SNR/10);%信噪比
% deltax = 0.1;
% x1=0:deltax:x0;
% x2=x0:deltax:N-1;
% step1=1500;
% step2=1000;
% 
% p_z_x_left=zeros(1,length(x1));
% p_z_x_right=zeros(1,length(x2));
% 
% N0=2*alpha^2/rho2;
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% for ii=1:length(x1)
%     p_z_x_left(1,ii)=exp(rho2/2*x1(ii)-step1)*...
%                      besseli(0,rho2*abs(N-x0+1/alpha*mean(wn.*(N-x1(ii))))-step2);
% end
% p_x_z_left=p_z_x_left/sum(p_z_x_left*deltax);
% 
% for jj=1:length(x2)
%     p_z_x_right(1,jj)=exp(rho2/2*x2(jj)-step1)*...
%                       besseli(0,rho2*abs(N-x2(jj)+1/alpha*mean(wn.*(N-x2(jj))))-step2);
%     if(p_z_x_right(1,jj)==inf)
%         number=randi(2)-1;
%         p_z_x_right(1,jj)=p_z_x_right(1,jj-1)*10^(6+number);
%     end
% end
% p_x_z_right=p_z_x_right/sum(p_z_x_right*deltax);
% plot(x1, p_x_z_left,'LineWidth',1);hold on
% plot(x2, p_x_z_right,'LineWidth',1);hold on
