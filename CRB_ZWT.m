clear all;clc;
%%
% 积分方式求CRB,不对

% figure;
% SNR = -15:2:15;
% rho2 = 10.^(SNR/10);%信噪比
% N =16;
% 
% deltax = 0.05;
% 
% rf = 10;
% e=2.718;
% 
% EEB=zeros(1,length(SNR));
% Num =100;
% SIGMAEE=zeros(1,length(SNR));
% crb=zeros(1,length(SNR));
% 
% x10 =N/2;
% x1 = 0: deltax : N-deltax;
% x2=0:deltax:x10;
% x3=x10:deltax:N;
% w=rf:1:200;    
% for snr=1:length(SNR)
%    N0=2*rf^2/rho2(snr); 
%    h=zeros(1,Num);
%    for run=1:Num
% %    p_z_x = zeros(1,length(x1));
%    p_z_x1 = zeros(1,length(x2));
%    p_z_x2 = zeros(1,length(x3));
%    wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
%      for kx=1:length(x2)             
%          A0=rho2(snr)-log(besseli(0,2*rho2(snr)*(1/rf)*abs(mean(wn))));
%          B0=(N-x10)*log(besseli(0,2*rho2(snr)*(1+(1/rf)*abs(mean(wn)))))+x10*log(besseli(0,2*rho2(snr)*(1/rf)*abs(mean(wn))));
%          p_z_x1(1,kx)=exp(A0*x1(kx)+B0);%x<x0
%      end
%      p_x_z1 = p_z_x1./sum(p_z_x1*deltax);
%      h1(1,run)=-sum(p_x_z1.*log(p_x_z1+eps))*deltax;
%      for jx=1:length(x3)
%          A1=rho2(snr)-log(besseli(0,2*rho2(snr)*(1+(1/rf)*abs(mean(wn)))));
%          B1=N*log(besseli(0,2*rho2(snr)*(1+(1/rf)*abs(mean(wn)))));
%          p_z_x2(1,jx)=exp(A1*x1(jx)+B1);%x>x0
%      end
%      p_x_z2 = p_z_x2./sum(p_z_x2*deltax);
%      h2(1,run)=-sum(p_x_z2.*log(p_x_z2+eps))*deltax;
%      h(1,run)=h1(1,run)+h2(1,run);
%      
%      A=0;
% for kx=rf:length(w)
%     A=A+log(besseli(0,2.*rf./N0.*w(kx))).*((w(kx)-rf)./N0.*exp(-((w(kx)-rf).^2)./(2.*N0))).*deltax;
% end
%     B=0;
% for kx=rf:length(w)
%     B=B+(log(besseli(0,2.*rf./N0.*w(kx))).^2).*((w(kx)-rf)./N0.*exp(-((w(kx)-rf).^2)./(2.*N0))).*deltax;
% end   
%    f=rho2(snr)-2.*rho2(snr).*A+B; 
%    
%    end
%    SIGMAEE(1,snr)=2^(2*mean(h))/(2*pi*e);
%    crb(1,snr)=1./f;
% end
% 
% semilogy(SNR,SIGMAEE,'LineWidth',1);hold on;
% semilogy(SNR,crb,'LineWidth',1);hold on;



%%
%熵误差和克拉美罗界
% clear all;clc;
% figure;
% SNR = -15:1:15;
% rho2 = 10.^(SNR/10);%信噪比
% N =16;
% 
% deltax = 0.05;
% 
% rf = 1;
% e=2.718;
% 
% EEB=zeros(1,length(SNR));
% Num =100;
% SIGMAEE=zeros(1,length(SNR));
% crb=zeros(1,length(SNR));
% 
% x10 =N/2;
% x1 = 0: deltax : N-deltax;
% x2=0:deltax:x10;
% x3=x10:deltax:N;
%     
% for snr=1:length(SNR)
%    N0=2*rf^2/rho2(snr); 
%    h=zeros(1,Num);
%    for run=1:Num
% %    p_z_x = zeros(1,length(x1));
%    p_z_x1 = zeros(1,length(x2));
%    p_z_x2 = zeros(1,length(x3));
%    wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
%      for kx=1:length(x2)             
%          A0=rho2(snr)-log(besseli(0,2*rho2(snr)*(1/rf)*abs(mean(wn))));
%          B0=(N-x10)*log(besseli(0,2*rho2(snr)*(1+(1/rf)*abs(mean(wn)))))+x10*log(besseli(0,2*rho2(snr)*(1/rf)*abs(mean(wn))));
%          p_z_x1(1,kx)=exp(A0*x1(kx)+B0);%x<x0
%      end
%      p_x_z1 = p_z_x1./sum(p_z_x1*deltax);
%      p_x_z1(isnan(p_x_z1))=0;
%      h1(1,run)=-sum(p_x_z1.*log(p_x_z1+eps))*deltax;
%      for jx=1:length(x3)
%          A1=rho2(snr)-log(besseli(0,2*rho2(snr)*(1+(1/rf)*abs(mean(wn)))));
%          B1=N*log(besseli(0,2*rho2(snr)*(1+(1/rf)*abs(mean(wn)))));
%          p_z_x2(1,jx)=exp(A1*x1(jx)+B1);%x>x0
%          m1 = max(p_z_x2(~isinf(p_z_x2)));
%          p_z_x2(find(p_z_x2==inf))=m1 ;
% 
%      end
%     
%      p_x_z2 = p_z_x2./sum(p_z_x2*deltax);
%      h2(1,run)=-sum(p_x_z2.*log(p_x_z2+eps))*deltax;
%      h(1,run)=h1(1,run)+h2(1,run);
%     
%      
%      E1 = zeros(1,length(x2));
%      E2 = zeros(1,length(x3));
%      a=(rho2(snr)-log(besseli(0,2.*rf./N0.*abs(rf+wn)))).^2;
%      E1(1,run)=mean(a);
%      b=(rho2(snr)-log(besseli(0,2.*rf./N0.*abs(wn)))).^2;
%      E2(1,run)=mean(b);
%      E(1,run)=E1(1,run)+E2(1,run);
%      
% 
%      
%    end
%    SIGMAEE(1,snr)=2^(2*mean(h))/(2*pi*e);
%    crb(1,snr)=1./sum(E);
% end
% 
% semilogy(SNR,SIGMAEE,'LineWidth',1);hold on;
% semilogy(SNR,crb,'--','LineWidth',1);hold on;
% xlabel('SNR');ylabel('SIGMAEE');
% legend('熵误差','克拉美罗界');
% title('恒模散射目标的克拉美罗界和熵误差')
%%
clear all;clc;
figure;
SNR = -10:1:15;
rho2 = 10.^(SNR/10);%信噪比
N =20;

deltax = 0.05;

rf = 1;
e=2.718;

EEB=zeros(1,length(SNR));
Num =100;
SIGMAEE=zeros(1,length(SNR));
crb=zeros(1,length(SNR));

x10 =N/2;
x1 = 0: deltax : N-deltax;
x2=0:1:x10/deltax;
x3=x10/deltax-1:1:N/deltax;
  
for snr=1:length(SNR)
   N0=2*rf^2/rho2(snr);
   
   h=zeros(1,Num);
   for run=1:Num
   p_z_x = zeros(1,length(x1));
   p_z_x1 = zeros(1,length(x2));
   p_z_x2 = zeros(1,length(x3));
   wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
     for kx=1:length(x2)             
         A0=rho2(snr)-log(besseli(0,2*rho2(snr)*(1/rf)*abs(mean(wn))));
         B0=(N-x10)*log(besseli(0,2*rho2(snr)*abs(1+mean(wn)./rf)))+x10*log(besseli(0,2*rho2(snr)*(1/rf)*abs(mean(wn))));
         p_z_x1(1,kx)=exp(A0*x1(kx)+B0);%x<x0
     end
     p_x_z1 = p_z_x1./sum(p_z_x1*deltax);
     p_x_z1(isnan(p_x_z1))=0;
     h1(1,run)=-sum(p_x_z1.*log(p_x_z1+eps))*deltax;
     for jx=1:length(x3)
         A1=rho2(snr)-log(besseli(0,2*rho2(snr)*abs(1+mean(wn)./rf)));
         B1=N*log(besseli(0,2*rho2(snr)*abs(1+mean(wn)./rf)));
         p_z_x2(1,jx)=exp(A1*x1(jx)+B1);%x>x0
         m1 = max(p_z_x2(~isinf(p_z_x2)));
         p_z_x2(find(p_z_x2==inf))=m1 ;

     end
     p_x_z2 = p_z_x2./sum(p_z_x2*deltax);
     h2(1,run)=-sum(p_x_z2.*log(p_x_z2+eps))*deltax;
     h(1,run)=h1(1,run)+h2(1,run);
     

     a=(rho2(snr)-log(besseli(0,2.*rf./N0.*abs(rf+wn)))).^2;
     b=(rho2(snr)-log(besseli(0,2.*rf./N0.*abs(wn)))).^2;
     E(1,run)=(mean(a)*N+mean(b)*N)/2;

     
    
     
     

     
   end
   SIGMAEE(1,snr)=2^(2*mean(h))/(2*pi*e);
   crb(1,snr)=1./mean(E);
   
  
   
   
end

semilogy(SNR,SIGMAEE,'LineWidth',1);hold on;
semilogy(SNR,crb,'--','LineWidth',1);hold on;
xlabel('SNR');ylabel('SIGMAEE');
legend('熵误差','克拉美罗界');


