% %%不同时间带宽积恒模散射目标的距离信息与SNR的关系
% clc;clear;close all;
% e=2.718;Pi=3.14159;
% % N=1024;
% for N=[1024,2048,4096]
% x0=N/2;n=1:N;
% SNR = -20:2:20;
% alpha=1;
% rho2 = 10.^(SNR/10);%信噪比
% deltax = 0.1;
% x1=0:deltax:N-deltax;
% x2=0:deltax:x0-deltax;
% x3=x0:deltax:N-deltax;
% Num =100;
% 
% 
% for snr=1:length(SNR)
%     N0=alpha^2/rho2(snr);
%      for run=1:Num
%         count1=0;
%         count2=0;
%         count3=0;
%         step=3.75;
%         L=1;
%         wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
%         a=rho2(snr)-log(besseli(0,2.*rho2(snr).*abs(mean(wn)./alpha)));
%         b=(x0+1).*log(besseli(0,2.*rho2(snr).*abs(mean(wn)./alpha)))+...
%           (N-x0).*log(besseli(0,2.*rho2(snr).*abs(1+mean(wn)./alpha)));
%         c=rho2(snr)-log(besseli(0,2.*rho2(snr).*abs(1+mean(wn)./alpha)));
%         d=N.*log(besseli(0,2.*rho2(snr).*abs(1+mean(wn)./alpha)));
%         %x<x0
%         for kx=1:length(x2)
%               p_z_x1(1,kx)=exp(a.*x2(kx)+b)/L;
%         end
%         p_x_z1 = p_z_x1./sum(p_z_x1*deltax); 
%         %x>x0
%         for kx=1:length(x3)
%               p_z_x2(1,kx)=exp(c*x3(kx)+d)/L;
%         end
%         p_x_z2 = p_z_x2./sum(p_z_x2*deltax); 
%         h(1,run)=-sum(p_x_z1.*log2(p_x_z1+eps)*deltax)...
%                  -sum(p_x_z2.*log2(p_x_z2+eps)*deltax); %理论熵，eps防止出现log0
%      end
%          I(1,snr)=log2(N)-mean(h);
% end
% plot(SNR,I,'LineWidth',2);
% hold on
% end
% legend('N=1024','N=2048','N=4096');
% xlabel('SNR');ylabel('I(Z;X)');
% title('不同时间带宽积恒模散射目标的距离信息与SNR的关系');


% %% 熵误差2.0
% clc;clear;close all;
% e=2.718;Pi=3.14159;
% N=16;
% x0=N/2;n=1:N;
% SNR = -20:1:20;
% alpha=1;
% rho2 = 10.^(SNR/10);%信噪比
% deltax = 0.01;
% x1=0:deltax:N-deltax;
% x2=0:deltax:x0-deltax;
% x3=x0:deltax:N-deltax;
% Num =100;
% 
% for i=1:length(SNR)
%     N0=2*alpha^2/rho2(i);
%      for run=1:Num
%         s1=0;
%         s2=0;
%         p_z_x1 = zeros(1,length(x2));
%         p_z_x2 = zeros(1,length(x3));
%         crb1=zeros(1,Num);
%         crb2=zeros(1,Num);
%         wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
%         a=rho2(i)-log(besseli(0,2.*rho2(i).*abs(mean(wn)./alpha)));
%         b=(x0).*log(besseli(0,2.*rho2(i).*abs(mean(wn)./alpha)))+...
%           (N-x0).*log(besseli(0,2.*rho2(i).*abs(1+mean(wn)./alpha)));
%         c=rho2(i)-log(besseli(0,2.*rho2(i).*abs(1+mean(wn)./alpha)));
%         d=N.*log(besseli(0,2.*rho2(i).*abs(1+mean(wn)./alpha)));
%         %x<x0
%         for kx=1:length(x2)
%               p_z_x1(1,kx)=exp(a.*x2(kx)+b);
%         end
%         p_x_z1 = p_z_x1./sum(p_z_x1*deltax); 
%         crb1(1,run)=(a).^2;
%         %x>x0
%         for kx=1:length(x3)
%               p_z_x2(1,kx)=exp(c*x3(kx)+d);
%         end
%         p_x_z2 = p_z_x2./sum(p_z_x2*deltax); 
%         crb2(1,run)=(c).^2;          
%         h(1,run)=-sum(p_x_z1.*log2(p_x_z1+eps)*deltax)...
%                  -sum(p_x_z2.*log2(p_x_z2+eps)*deltax); %理论熵，eps防止出现log0
%         crb(1,run)=crb1(1,run);
%         crb(1,Num+run)=crb2(1,run);
%      end
%      CRB(1,i)=1./mean(crb);
%      SIGMAEE(1,i)=2^(2*mean(h))/(2*pi*e);
% end
% semilogy(SNR,SIGMAEE,'LineWidth',2);
% hold on;grid on;
% semilogy(SNR,CRB,'LineWidth',2);
% hold on
% xlabel('SNR');ylabel('ERROR');
% legend('EE','CRB');
% title('熵误差');


%% EE CRB
clc;clear;close all;
e=2.718;
N=16;
x0=N/2;n=1:N;
SNR = -10:1:25;
alpha=1;
rho2 = 10.^(SNR/10);%信噪比
deltax = 0.01;
x1=0:deltax:N-deltax;
x2=0:deltax:x0-deltax;
x3=x0:deltax:N-deltax;
Num =100;
step1=1e1;
step2=5;
for i=1:length(SNR)
    snr=rho2(i);
    N0=2*alpha^2/snr;
     for run=1:Num
        s1=0;
        s2=0;
        p_z_x1 = zeros(1,length(x2));
        p_z_x2 = zeros(1,length(x3));
        crb1=zeros(1,Num);
        crb2=zeros(1,Num);
        wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
        a=(snr-log(besseli(0,2.*snr.*abs(mean(wn)./alpha))));
        b=((x0).*log(besseli(0,2.*snr.*abs(mean(wn)./alpha)))+...
          (N-x0).*log(besseli(0,2.*snr.*abs(1+mean(wn)./alpha))))/step1;
        c=(snr-log(besseli(0,2.*snr.*abs(1+mean(wn)./alpha))));
        d=N.*log(besseli(0,2.*snr.*abs(1+mean(wn)./alpha)))/step1;
        %x<x0
        for kx=1:length(x2)
              p_z_x1(1,kx)=exp((a.*x2(kx)+b)/step2);
%               m1 = max(p_z_x1(~isinf(p_z_x1)));
%               p_z_x1(find(p_z_x1==inf))=m1 ;
        end
        p_x_z1 = p_z_x1./sum(p_z_x1*deltax);
%         p_x_z1(isnan(p_x_z1))=0;
        crb1(1,run)=(a).^2;
        %x>x0
        for kx=1:length(x3)
              p_z_x2(1,kx)=exp((c*x3(kx)+d)/step2);
%               m2 = max(p_z_x2(~isinf(p_z_x2)));
%               p_z_x2(find(p_z_x2==inf))=m2 ;
        end
        p_x_z2 = p_z_x2./sum(p_z_x2*deltax); 
        crb2(1,run)=(c).^2;          
        h(1,run)=-sum(p_x_z1.*log2(p_x_z1+eps)*deltax)...
                 -sum(p_x_z2.*log2(p_x_z2+eps)*deltax); %理论熵，eps防止出现log0
%         crb(1,run)=crb1(1,run);
%         crb(1,Num+run)=crb2(1,run);
        crb(1,run)=(N*crb1(1,run)+N*crb1(1,run))/2;
     end
     I(1,i)=log2(N)-mean(h);
     CRB(1,i)=1./mean(crb);
     EE(1,i)=2^(2*mean(h))/(2*pi*e);
end
plot(SNR,I,'LineWidth',2);hold on;grid on;
% semilogy(SNR,CRB,'--','LineWidth',2);hold on
% semilogy(SNR,EE,'LineWidth',2);hold on;grid on;
% semilogy(SNR,CRB,'--','LineWidth',2);hold on
% xlabel('SNR');ylabel('ERROR');
% legend('EE','CRB');




% %% 熵误差1.0
% clc;clear;close all;
% e=2.718;Pi=3.14159;
% N=16;
% x0=N/2;n=1:N;
% SNR = -20:1:20;
% alpha=1;
% rho2 = 10.^(SNR/10);%信噪比
% deltax = 0.1;
% x1=0:deltax:N-deltax;
% x2=0:deltax:x0-deltax;
% x3=x0:deltax:N-deltax;
% Num =100;
% 
% for snr=1:length(SNR)
%     N0=alpha^2/rho2(snr);
%      for run=1:Num
%         count1=0;
%         count2=0;
%         count3=0;
%         step=3.75;
%         L=1;
%         wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
%         a=rho2(snr)-log(besseli(0,2.*rho2(snr).*abs(mean(wn)./alpha)));
%         b=(x0+1).*log(besseli(0,2.*rho2(snr).*abs(mean(wn)./alpha)))+...
%           (N-x0).*log(besseli(0,2.*rho2(snr).*abs(1+mean(wn)./alpha)));
%         c=rho2(snr)-log(besseli(0,2.*rho2(snr).*abs(1+mean(wn)./alpha)));
%         d=N.*log(besseli(0,2.*rho2(snr).*abs(1+mean(wn)./alpha)));
%         %x<x0
%         for kx=1:length(x2)
%               p_z_x1(1,kx)=exp(a.*x2(kx))/L;
% %                   if p_z_x1(1,kx)>1e300
% %                       count1=count1+step;
% %                       p_z_x1(1,kx)=randi(9)*10^7*10^count1;
% %                   end
%         end
%         p_x_z1 = p_z_x1./sum(p_z_x1*deltax); 
%         crb1(1,run)=(a.^2);
%         %x>x0
%         for kx=1:length(x3)
%               p_z_x2(1,kx)=exp(c*x3(kx))/L;
% %                   if p_z_x2(1,kx)>1e300
% %                       count2=count2+step;
% %                       p_z_x2(1,kx)=randi(9)*10^10/10^count2;
% %                   end
%         end
%         p_x_z2 = p_z_x2./sum(p_z_x2*deltax); 
%         crb2(1,run)=(c.^2);          
%         h(1,run)=-sum(p_x_z1.*log2(p_x_z1+eps)*deltax)...
%                  -sum(p_x_z2.*log2(p_x_z2+eps)*deltax); %理论熵，eps防止出现log0
%         crb(1,run)=crb1(1,run);
%         crb(1,Num+run)=crb2(1,run);
%      end
%      CRB(1,snr)=1./mean(crb);
%      SIGMAEE(1,snr)=2^(2*mean(h))/(2*pi*e);
% %      if SIGMAEE(1,snr)==CRB(1,snr)
% %         SIGMAEE(1,snr:length(SNR))=CRB(1,snr:length(SNR));
% %      end
% 
% %          I(1,snr)=log2(N)-mean(h);
% end
% %     plot(SNR,I,'LineWidth',2);
% %     hold on
% semilogy(SNR,SIGMAEE,'LineWidth',2);
% hold on;grid on;
% semilogy(SNR,CRB,'--','LineWidth',2);
% hold on
% xlabel('SNR');ylabel('ERROR');
% legend('EE','CRB');
% title('熵误差');






% %% 微分熵 熵误差 integral函数积分
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;e=2.718;Pi=3.14159;
% N=16;x0=N/2;n=1:N;N0=k*T;
% 
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% w=mean(wn);
% for SNR=0.01:0.01:10
%     alpha=sqrt(SNR*N0);
%     a=SNR-log(besseli(0,2.*SNR.*abs(w./alpha)));
%     b=(x0+1).*log(besseli(0,2.*SNR.*abs(w./alpha)))+...
%       (N-x0).*log(besseli(0,2.*SNR.*abs(1+w./alpha)));
%     c=SNR-log(besseli(0,2.*SNR.*abs(1+w./alpha)));
%     d=N.*log(besseli(0,2.*SNR.*abs(1+w./alpha)));
%     
%     in=@(x) exp(SNR.*x).*(besseli(0,2*SNR*abs(w/a))).^(x0-x+1).*(besseli(0,2*SNR*abs(1+w/a))).^(N-x0+1);
%     A1=integral(in,0,x0);%x<x0的分母
%     in=@(x) exp(SNR.*x).*(besseli(0,2*SNR*abs(1+w/a))).^(N-x+1);
%     A2=integral(in,x0,N);%x>x0的分母
%     
%     in1=@(x) (a*x+b).*exp(a*x+b)/(-A1)+exp(a*x+b)*log(A1)/(A1);
%     s1=integral(in1,0,x0);
%     in2=@(x) (c*x+d).*exp(c*x+d)/(-A2)+exp(c*x+d)*log(A2)/(A2);
%     s2=integral(in2,x0,N-1);
%     s=s1+s2;
%     figure(1)
%     semilogy(SNR,s,'r.','LineWidth',2);
%     hold on;grid on; 
%     figure(2)
%     y=2^(2*s)/(2*Pi*e);
%     semilogy(SNR,y,'r.','LineWidth',4);
%     hold on
% end
% ylabel('微分熵');
% xlabel('SNR/dB');

% % set(gca,'Ytick',[1.0e-7 1.0e-6 1.0e-5 1.0e-4 1.0e-3 1.0e-2 1.0e-1 1.0e0 1.0e1 1.0e10]);
% % % set(gca,'Ytick',[0.0001 0.001 0.01 0.1 1.0 10.0]);
% % set(gca,'YTickLabel',{'10','1','0.1','0.01','0.001','0.0001'});
% % ylim([0 15]);
% % y_tick={'1.00e-06','1.00e-05','1.00e-04','1.00e-03','1.00e-02','1.00e-01','1.00e-00','1.00e1'};
% % set(gca,'YtickLabel',y_tick)
% % figure(1)
% % set(gca,'Ytick',[1.0e-7 1.0e-6 1.0e-5 1.0e-4 1.0e-3 1.0e-2 1.0e-1 1.0e0 1.0e1 1.0e10]);
% % ylabel('微分熵');
% % xlabel('SNR/dB');
% % figure(2)
% % set(gca,'Ytick',[1.0e-7 1.0e-6 1.0e-5 1.0e-4 1.0e-3 1.0e-2 1.0e-1 1.0e0 1.0e1 1.0e10]);
