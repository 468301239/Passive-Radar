%%除下来修改后2.0
clc;clear;close all;
k=1.380649*10^-23;T=273.15+25;deltax = 0.01;
N0=k*T;
pi=3.14;%SSNR=a/N0*[0.001 0.01 0.1];
p1=0:0.01:1;p0=1-p1;
SNR=1;
a=sqrt(SNR*N0);
%n=1:NN;
% wn=sqrt(N0/2).*(randn(1)+1i*randn(1));
% wn=sqrt(N0/2).*randn(N,1);
% n=1:NN;
% wn=zeros(NN,NN);
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n))); 
in=1;Tao=0;
for N=[256 512 1024]
    n=1:N;
    wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n))); 
    for x=1:N
        for kx=x:N
            w=wn(kx);
            in=in*((besseli(0,2*a/N0.*abs(w)))/(exp(SNR)));
        end
%             Tao=sum(in*deltax);
        Tao=Tao+in;
        pfa=(p1*Tao)./(N*p0+p1*Tao);
    end
    plot(p1,pfa,'--','LineWidth',1);
    hold on;
end
a=0:10^-4:1;
b=a;
plot(a,b,'r','LineWidth',1);
xlabel('P(1)');ylabel('P_F_A');
legend( 'N=256','N=512','N=1024','P_F_A=P(1)');
set(gca,'FontName','Times New Roman','FontSize',12)



% %%
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;
% N=128;N0=100*k*T;a=N0;pi=3.14;%SSNR=a/N0*[0.001 0.01 0.1];
% SSNR=a/N0*[0.01];
% SNR=0.01*a/N0;n=1:N;
% wn=sqrt(N0/2).*(randn(1)+1i*randn(1));
% % wn=sqrt(N0/2).*randn(N,1);
% for i=1:length(SSNR)
%     SNR=SSNR(i);
%     T=@(x) exp(-(N-x+1).*SNR).*(besseli(0,2*a/N0.*abs(wn))).^(N-x+1);
%     Tao=integral(T,0,N)/N;
% 
%     p1=0:0.01:1;
%     pfa=(p1*Tao)./(1-p1+p1*Tao);
%     plot(p1,pfa,'--','LineWidth',1);
%     hold on;
% end
% a=0:10^-4:1;
% b=a;
% plot(a,b,'r','LineWidth',1);
% xlabel('P(1)');ylabel('P_F_A');
% legend( 'N=128','P_F_A=P(1)');
% set(gca,'FontName','Times New Roman','FontSize',12)

% function[wn]=Wn(x)
% k=1.380649*10^-23;T=273.15+25;
% N0=100*k*T;
% wn=sqrt(N0/2)*(rand(fix(x))+1i*rand(fix(x)));
% end


% %%
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;
% NN=128;N0=100*k*T;a=N0;pi=3.14;%SSNR=a/N0*[0.001 0.01 0.1];
% SNR=0.01*a/N0;n=1:NN;
% p1=0:10^-4:1; 
% for N=NN/4:NN/4:NN
%     for run=1:N
%         wn(:,run)=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n))); 
%         w=sum(wn(:,run));
%         T=@(x) exp(-(N-x+1).*SNR).*(besseli(0,2*a/N0.*abs(w))).^(N-x+1);
%         Tao=integral(T,1,N)/N;
%         pfa=(p1*Tao)./(1-p1+p1*Tao);
%         plot(p1,pfa,'--','LineWidth',1);
%         hold on;
%     end
% end
% % legend( 'N=32', 'N=64','N=96','N=128');
% a=0:10^-4:1;
% b=a;
% plot(a,b,'r-','LineWidth',1);
% xlabel('P(1)');ylabel('P_F_A');
% legend( 'N=32', 'N=64','N=96','N=128','P_F_A=P(1)');
% set(gca,'FontName','Times New Roman','FontSize',12)


% %% 除下来修改后
% % clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;
% N=1024;N0=100*k*T;a=N0;pi=3.14;%SSNR=a/N0*[0.001 0.01 0.1];
% SNR=0.001*a/N0;n=1:N;
% p1=0:10^-4:1; 
% for run=1:N
%     wn(:,run)=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n))); 
%     w=sum(wn(:,run));
% %     T=@(x) exp(-(N-x+1).*SNR).*(besseli(0,2*a/N0.*abs(w))).^(N-x+1);
%     T=@(x) ((besseli(0,2*a/N0.*abs(w))).^(N-x+1)/(exp(-SNR))).^(N-x);
%     Tao=integral(T,1,N)/N;
%     pfa=(p1*Tao)./(1-p1+p1*Tao);
%     plot(p1,pfa,'--','LineWidth',1);
%     hold on;
% end
% a=0:10^-4:1;
% b=a;
% plot(a,b,'r','LineWidth',1);
% xlabel('P(1)');ylabel('P_F_A');
% % legend( 'N=32', 'N=64','N=96','N=128','P_F_A=P(1)');
% set(gca,'FontName','Times New Roman','FontSize',12)

% %%除下来修改后
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;
% NN=1024;N0=100*k*T;a=N0;pi=3.14;%SSNR=a/N0*[0.001 0.01 0.1];
% SNR=0.01*a/N0;n=1:NN;
% wn=sqrt(N0/2).*(randn(1)+1i*randn(1));
% % wn=sqrt(N0/2).*randn(N,1);
% for N=NN/4:NN/4:NN
%     T=@(x) ((besseli(0,2*a/N0.*abs(wn))).^(N-x)/(exp(-SNR))).^(N-x);
%     Tao=integral(T,0,N-1)/N;
% 
%     p1=0:0.01:1;
%     pfa=(p1*Tao)./(1-p1+p1*Tao);
%     plot(p1,pfa,'--','LineWidth',1);
%     hold on;
% end
% a=0:10^-4:1;
% b=a;
% plot(a,b,'r','LineWidth',1);
% xlabel('P(1)');ylabel('P_F_A');
% legend( 'N=256','N=512','N=768','N=1024','P_F_A=P(1)');
% set(gca,'FontName','Times New Roman','FontSize',12)

