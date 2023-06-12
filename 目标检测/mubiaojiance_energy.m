%% IT NP检测信息对比
clc;clear;close all;
N=500;
deltax=0.05;                     
x=0:deltax:N-deltax;
p1=0:0.0001:1;
p0=1-p1;
rho2=[0];
size=length(rho2);
for i=1:size
    T_1_z1_up=0;
    T_1_z1_down=0;
    T_1_z1=0;   
    T_1_z0=0;
    pd_NP=0;
    SNR=10^(rho2(i)/10);    
    for kx=1:length(x)
        T_1_z1_up=T_1_z1_up+1/sqrt(N+2*(N-x(kx))*SNR)*deltax;
        T_1_z1_down=T_1_z1_down+1/sqrt(N)*exp(-(N-x(kx))^2*SNR^2/2/N)*deltax;
%         T_1_z1=T_1_z1+1/sqrt(N+2*(N-x(kx))*SNR)./(1/sqrt(N)*exp(-(N-x(kx))^2*SNR^2/2/N))*deltax;
        T_1_z0=T_1_z0+1/N*sqrt(N/(N+2*(N-x(kx))*SNR))*exp(-(N-x(kx))^2*SNR^2./(2*N+4*(N-x(kx))*SNR))*deltax;
        pfa_NP=p1.*T_1_z0./(p0+p1.*T_1_z0);
        pd_NP=pd_NP+0.5*erfc((sqrt(2*N)*erfcinv(2*pfa_NP)-SNR*(N-x(kx)))/(sqrt(2*N+4*(N-x(kx))*SNR)))/N*deltax; 
    end
    T_1_z1=T_1_z1_up./T_1_z1_down;
    p1z0=p1.*T_1_z0./(p0+p1.*T_1_z0);
    p1z1=p1.*T_1_z1./(p0+p1.*T_1_z1);
    p0z0=p0./(p0+p1.*T_1_z0);
%     p0z0=p0./(p0+p1.*T_1_z1);
%     p1z1=ones(1,10001);
%     p0z0=ones(1,10001);
    p0z1=p0./(p0+p1.*T_1_z1);
    I_IT=-p1.*log2(p1)-p0.*log2(p0)...
         +p1z0.*log2(p1z0)...
         +p1z1.*log2(p1z1)...
         +p0z0.*log2(p0z0)...
         +p0z1.*log2(p0z1);
     A=p1.*pd_NP./(p1.*pd_NP-p1.*pfa_NP+pfa_NP);
     D=(1-p1).*(1-pfa_NP)./(1-p1.*pd_NP+p1.*pfa_NP-pfa_NP);
     I_NP=-p1.*log2(p1)-p0.*log2(p0)...
          -(p1.*pd_NP-p1.*pfa_NP+pfa_NP).*(-A.*log2(A)-(1-A).*log2(1-A))...
          -(1-p1.*pd_NP+p1.*pfa_NP-pfa_NP).*(-D.*log2(D)-(1-D).*log2(1-D));
    h=figure;
    semilogx(p1,I_IT,'r',p1,I_NP,'--','LineWidth',2);
    hold on;grid on;
%     semilogx(p1,I_NP,'--','LineWidth',2);
end
set(gca,'FontSize',12)
xlabel('p1');ylabel('I');
% title('IT NP 检测信息对比  IT检测概率统计量积分');
% legend('SNR=   0dB  IT','SNR=   0dB  NP');
saveas(h,'IT NP检测信息对比','fig')


% %% IT NP PFA PD 先积去x ROC对比 一个图
% clc;clear;close all;
% N=64;
% deltax=0.05;                     
% x=0:deltax:N-deltax;
% pfa=0:0.00001:1;
% % Pi=3.1415926;
% rho2=[0];
% size=length(rho2);
% IT_NP_ROC=figure;
% for i=1:size
%     T_1_z1_up=0;
%     T_1_z1_down=0;
%     T_1_z0=0;
%     pd_NP=0;
%     SNR=10^(rho2(i)/10);    
% %     alpha2=sqrt(SNR*N0);
%     for kx=1:length(x)
%         T_1_z1_up=T_1_z1_up+1/sqrt(N+2*(N-x(kx))*SNR)*deltax;
%         T_1_z1_down=T_1_z1_down+1/sqrt(N)*exp(-(N-x(kx))^2*SNR^2/2/N)*deltax;
%         T_1_z0=T_1_z0+1/N*sqrt(N/(N+2*(N-x(kx))*SNR))*exp(-(N-x(kx))^2*SNR^2/(2*N+4*(N-x(kx))*SNR))*deltax;
%         pd_NP=pd_NP+0.5*erfc((sqrt(2*N)*erfcinv(2*pfa)-SNR*(N-x(kx)))/(sqrt(2*N+4*(N-x(kx))*SNR)))/N*deltax; 
%     end
%     T_1_z1=T_1_z1_up./T_1_z1_down;
%     pd_IT=pfa.*T_1_z1./((1-pfa).*T_1_z0+pfa.*T_1_z1);
%     semilogx(pfa,pd_IT,pfa,pd_NP,'--','LineWidth',2);
%     hold on;grid on;
% end
% set(gca,'FontSize',12)
% xlabel('P_F_A');ylabel('P_D');
% title('IT NP PFA PD 先积去x ROC对比');
% % legend('SNR=-20dB  IT','SNR=-20dB  NP',...
% %        'SNR= -10dB  IT','SNR= -10dB  NP',...
% %        'SNR=  -5dB  IT','SNR=  -5dB  NP');
% saveas(IT_NP_ROC,'IT NP','fig')

% %% IT  NP PFA PD x没积去
% clc;clear;close all;
% N=512;                  
% x=0;
% N0=1;
% rho2=[ -10 -5 0];
% pfa=0:0.001:1;
% for i=1:length(rho2)
%     SNR=10^(rho2(i)/10);
%     alpha=sqrt(SNR*N0);
%     fx=exp((N+(N-x)*SNR)*((N-x)^2*SNR^2)/(N*(N+2*(N-x)*SNR)));
%     pd_IT=pfa.*fx./((fx-1).*pfa+1);
%     semilogx(pfa,pd_IT,'LineWidth',2);
%     hold on;
%     Th=sqrt(2*N*N0^2).*erfcinv(2*pfa)+N*N0;
%     pd_NP=0.5*erfc((sqrt(2*N)*erfcinv(2*pfa)-SNR*(N-x))/(sqrt(2*N+4*(N-x)*SNR))); 
%     semilogx(pfa,pd_NP,'--','LineWidth',2);
%     hold on;
% end
% set(gca,'FontName','Times New Roman','FontSize',12)
% legend('SNR=-10dB  IT','SNR=-10dB  NP',...
%        'SNR=  -5dB  IT','SNR=  -5dB  NP',...
%        'SNR=   0dB  IT','SNR=   0dB  NP');

% % IT NP PFA PD 先积去x ROC对比 多个图
% clc;clear;close all;
% N=64;
% deltax=0.01;                     
% x=0:deltax:N-1;
% pfa=0:0.001:1;
% Pi=3.1415926;
% N0=1;
% rho2=[-20 -10 -5 0];
% % rho2=[0 5 10 20];
% size=length(rho2);
% for i=1:size
%     subplot(ceil(size/2),fix(size/2),i);
%     T_1_z1_up=0;
%     T_1_z1_down=0;
%     T_1_z0=0;
%     pd_NP=0;
%     SNR=10^(rho2(i)/10);    
%     alpha2=sqrt(SNR*N0);
%     for kx=1:length(x)
%         T_1_z1_up=T_1_z1_up+1/sqrt(N+2*(N-x(kx))*SNR)*deltax;
%         T_1_z1_down=T_1_z1_down+1/sqrt(N)*exp(-(N-x(kx))^2*SNR^2/2/N)*deltax;
%         T_1_z0=T_1_z0+1/N*sqrt(N/(N+2*(N-x(kx))*SNR))*exp(-(N-x(kx))^2*SNR^2/(2*N+4*(N-x(kx))*SNR))*deltax;
%         pd_NP=pd_NP+0.5*erfc((sqrt(2*N)*erfcinv(2*pfa)-SNR*(N-x(kx)))/(sqrt(2*N+4*(N-x(kx))*SNR)))/N*deltax; 
%     end
%     T_1_z1=T_1_z1_up/T_1_z1_down;
%     pd_IT=pfa.*T_1_z1./(-T_1_z0*pfa+T_1_z0+pfa*T_1_z1);
%     semilogx(pfa,pd_IT,'r','LineWidth',1);
%     hold on;
%     semilogx(pfa,pd_NP,'b--','LineWidth',1);
%     hold on;grid on;
%     set(gca,'FontName','Times New Roman','FontSize',12)
%     title(sprintf('SNR=%d dB',rho2(i)));  
% end
% xlabel('P_F_A');ylabel('P_D');
% legend('IT','NP');


% %% IT PFA PD 先积去x
% clc;clear;close all;
% N=64;
% deltax=0.01;                     
% x=0:deltax:N-1;
% pfa=0:0.001:1;
% % Pi=3.1415926;
% N0=1;
% rho2=[-30 -20 -10 -5 0];
% for i=1:length(rho2)
%     T_1_z1_up=0;
%     T_1_z1_down=0;
%     T_1_z0=0;
%     SNR=10^(rho2(i)/10);    
% %     alpha2=sqrt(SNR*N0);
%     for kx=1:length(x)
%         T_1_z1_up=T_1_z1_up+1/sqrt(N+2*(N-x(kx))*SNR)*deltax;
%         T_1_z1_down=T_1_z1_down+1/sqrt(N)*exp(-(N-x(kx))^2*SNR^2/2/N)*deltax;
%         T_1_z0=T_1_z0+1/N*sqrt(N/(N+2*(N-x(kx))*SNR))*exp(-(N-x(kx))^2*SNR^2/(2*N+4*(N-x(kx))*SNR))*deltax;
%     end
%     T_1_z1=T_1_z1_up/T_1_z1_down;
%     pd=pfa.*T_1_z1./(-T_1_z0*pfa+T_1_z0+pfa*T_1_z1);
%     semilogx(pfa,pd,'LineWidth',2);
%     hold on;grid on;
% end
% set(gca,'FontSize',12)
% xlabel('P_F_A');ylabel('P_D');
% title('IT PFA PD 先积去x');
% legend('SNR=-30dB  IT','SNR=-20dB  IT','SNR=-10dB  IT','SNR= -5dB   IT','SNR=  0dB   IT');

% %% NP PFA PD 积去x，先积后积相同
% clc;clear;close all;
% N=64;
% deltax=0.05;                     
% x=0:deltax:N;
% pfa=0:0.001:1;
% rho2=[-30 -20 -10 -5 0];
% for i=1:length(rho2)
%     SNR=10^(rho2(i)/10);
%     pd=0;
%     for kx=1:length(x)
%         pd=pd+0.5*erfc((sqrt(2*N)*erfcinv(2*pfa)-SNR*(N-x(kx)))/(sqrt(2*N+4*(N-x(kx))*SNR)))/N*deltax; 
%     end
%     semilogx(pfa,pd,'LineWidth',2);
%     hold on;grid on;
% end
% set(gca,'FontSize',12)
% xlabel('P_F_A');ylabel('P_D');
% xlim([0 1]);
% title('NP PFA PD 先积去x');
% legend('SNR=-30dB  NP','SNR=-20dB  NP','SNR=-10dB  NP','SNR= -5dB   NP','SNR=  0dB   NP');




% %% IT PFA PD 后积去x
% clc;clear;close all;
% N=64;
% deltax=0.01;                     
% x=0:deltax:N;
% pfa=0:0.001:1;
% rho2=[-30 -20 -10 -5 0];
% fx=0;
% for i=1:length(rho2)
%     SNR=10^(rho2(i)/10);  
%     for kx=1:length(x)
%         fx=fx+exp((N+(N-x(kx))*SNR)*((N-x(kx))^2*SNR^2)/(N*(N+2*(N-x(kx))*SNR)))/N*deltax;
%     end
%     pd=pfa.*fx./((fx-1).*pfa+1);
%     semilogx(pfa,pd,'LineWidth',2);
%     hold on;grid on;
% end
% set(gca,'FontSize',12)
% xlabel('P_F_A');ylabel('P_D');
% title('IT PFA PD 后积去x');
% legend('SNR=-30dB  IT','SNR=-20dB  IT','SNR=-10dB  IT','SNR= -5dB   IT','SNR=  0dB   IT');

% %% IT NP 后积去时延x的ROC对比
% clc;clear;close all;
% N=64;
% deltax=0.01;                     
% x=0:deltax:N-1;
% pfa=0:deltax:1;
% i=1;                                
% rho2=[-5 0 10];
% for n=1:length(rho2)
%     SNR=10^(rho2(i)/10);
%     fx=0;pd_NP=0;
%     subplot(2,2,i);
%     for kx=1:length(x)
%         fx=fx+exp((N+(N-x(kx))*SNR)*((N-x(kx))^2*SNR^2)/(N*(N+2*(N-x(kx))*SNR)))/N*deltax;
%         pd_NP=pd_NP+0.5*erfc((sqrt(2*N)*erfcinv(2*pfa)-SNR*(N-x(kx)))/(sqrt(2*N+4*(N-x(kx))*SNR)))/N*deltax; 
%     end
%     pd_IT=pfa.*fx./((fx-1).*pfa+1);
%     semilogx(pfa,pd_IT,'LineWidth',1);
%     hold on;grid on;
%     semilogx(pfa,pd_NP,'--','LineWidth',1);
%     hold on;
%     set(gca,'FontName','Times New Roman','FontSize',12)
%     xlabel('P_F_A');ylabel('P_D');
%     title(sprintf('SNR=%d dB',rho2(i)));    
%     i=i+1;
% end
% legend('IT','NP');


% %% IT PFA PD x没积去
% clc;clear;close all;
% N=64;
% SNR=0.1;
% deltax=0.05;                     
% x=0;
% rho2=[-30 -20 -10 -5 0];
% pfa=0:0.001:1;
% for i=1:length(rho2)
%     SNR=10^(rho2(i)/10);  
%     fx=exp((N+(N-x)*SNR)*((N-x)^2*SNR^2)/(N*(N+2*(N-x)*SNR)));
%     pd=pfa.*fx./((fx-1).*pfa+1);
%     semilogx(pfa,pd,'LineWidth',2);
%     hold on;
% end
% set(gca,'FontName','Times New Roman','FontSize',12)
% legend('SNR=-30dB  IT','SNR=-20dB  IT','SNR=-10dB  IT','SNR= -5dB   IT','SNR=  0dB   IT');

