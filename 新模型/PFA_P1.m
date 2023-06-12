% close all;
% clc;tic 

a=0:10^-4:1;
b=a;
plot(a,b,'r')
hold on
grid on
%%
% d =0.01; %-20dB
% N = 2^7;
% P1=0:10^-4:1;
% P0=1-P1;
% N0=1;
% alpha=sqrt(N0*d);
% p=sqrt(1/(1+2*d))*exp(-N*d^2/(2*(1+2*d)));
% Pf0=P1*p./(P0+P1*p);
% plot(P1,Pf0
% hold on;
% grid on
%%
SNR=0.001;
rho_2 =SNR;
N = 2^4;
n = -N/2 : 1 : N/2-1;
deltax = 0.01;
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
    p(kx)=exp(-rho_2.*(N-x(kx))).*besseli(0,2*alpha/N0*abs(sum(wn(:,run))))^(N-x(kx));
     end
     P=sum(p*deltax);
    Pf0(run,:)=P1*P./(N*P0+P1*P);
end
Pf=mean(Pf0);
plot(P1,Pf)
hold on;
grid on

xlabel('\pi(1)');
ylabel('P_F_A');
legend('P_F_A=\pi(1)', 'N=16', 'N=32','N=64');
set(gca,'FontName','Times New Roman','FontSize',12)