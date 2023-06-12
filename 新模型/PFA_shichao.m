close all;clear;
clc;tic 

a=0:10^-4:1;
b=a;
plot(a,b,'r')
hold on
grid on

% SNR=0;
% rho_2 = 10.^(SNR/10);
% N = 2^7;
% n = -N/2 : 1 : N/2-1;
% deltax = 0.01;
% x = -N/2 : deltax : N/2-deltax;
% P1=0:10^-4:1;
% P0=1-P1;
% N0=1;
% num=100;
% alpha=sqrt(N0*rho_2);
% wn=zeros(N,num);
% Pf0=zeros(num,length(P1));
% p=zeros(1,length(x));
% for run=1:num
%     wn(:,run) = sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
%     for kx = 1 : length(x)
%     p(kx)=besseli(0,2*alpha/N0*abs(sum(wn(:,run)'*sinc(n'-x(kx)))));
%     end
%     P=exp(-rho_2)*sum(p*deltax);
%     Pf0(run,:)=P1*P./(N*P0+P1*P);
% end
% Pf=mean(Pf0);
% plot(P1,Pf,'kp','MarkerIndices',1:800:length(P1))
% hold on;
% grid on

SNR=10;
rho_2=SNR;
%  rho_2 = 10.^(SNR/10);
N = 2^8;
n = -N/2 : 1 : N/2-1;
deltax = 0.1;
x = -N/2 : deltax : N/2-deltax;
P1=0:10^-4:1;
P0=1-P1;
N0=1;
num=100;
alpha=sqrt(N0*rho_2);
wn=zeros(N,num);
Pf0=zeros(num,length(P1));
p=zeros(1,length(x));
for run=1:num
    wn(:,run) = sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
    for kx = 1 : length(x)
    p(kx)=besseli(0,2*alpha/N0*abs(sum(wn(:,run)'*sinc(n'-x(kx)))));
    end
    P=exp(-rho_2)*sum(p*deltax);
    Pf0(run,:)=P1*P./(N*P0+P1*P);
end
Pf=mean(Pf0);
plot(P1,Pf,'c')
hold on;
toc

xlabel('\pi(1)');
ylabel('P_F_A');
legend('P_F_A=\pi(1)', 'N=1024', 'N=2048');
set(gca,'FontName','Times New Roman','FontSize',12)