%% 检测信息
%% 恒模
% IP
% d =1;%信噪比
% N = 32;%时间带宽积
% pfa=0:0.0001:1;
% n = 0 : 1 : N-1;
% deltax = 1;
% x = 0 : deltax : N-deltax;
% % IP PFA PD
% for kx = 1 : length(x)
%     y(kx)=sqrt(1./(2*d*(N-x(kx))+N)).*exp(-(d*(N-x(kx)))^2./(2*(N+2*d*(N-x(kx)))));
%     v(kx)=1/sqrt(N+2*d*(N-x(kx)));
%     u(kx)=exp(-(d*(N-x(kx)))^2./(2*N));
% end
% y=sum(y);
% v=sum(v);
% u=sum(u);
% p1=(pfa.*sqrt(N))./(y+sqrt(N).*pfa-y.*pfa);
% p0=1-p1;
% PD=(p1.*v*sqrt(N))./(p1.*v*sqrt(N)+(1-p1).*u);
% %NP PFA PD
% z=zeros(1,10001);
% for kx = 1 : length(x)
%     z=z+1/2*erfc((sqrt(2*N)*erfcinv(2*p1)-d*(N-x(kx)))/(sqrt(2*(2*(N-x(kx))*d+N))));
% end
% PD1=z./N;
% 
% HV=-p0.*log(p0)-p1.*log(p1);
% 
% p1_z0=pfa;%PFA
% p0_z0 = 1 - p1_z0;
% hv_z0=-(p1_z0 .* log(p1_z0))-(p0_z0 .* log(p0_z0));
% 
% p1_z1=PD;%PD
% p0_z1 = 1 - p1_z1;  
% hv_z1 = -(p1_z1 .* log(p1_z1))-(p0_z1 .* log(p0_z1)) ;
% 
% IIT=HV-hv_z0-hv_z1;
% semilogx(p1,IIT,'g');%IT
% set(gca,'XScale','log')
% hold on;
% 
% %NP
% A=p1.*PD1./(p1.*PD1-p1.*pfa+pfa);
% D=(1-p1).*(1-pfa)./(1-p1.*PD1+p1.*pfa-pfa);
% INP=HV-(p1.*PD1-p1.*pfa+pfa).*(-A.*log(A)-(1-A).*log(1-A))-(1-p1.*PD1+p1.*pfa-pfa).*(-D.*log(D)-(1-D).*log(1-D));
% semilogx(p1,INP,'r');%np
% set(gca,'XScale','log')
% xlabel('PFA'); ylabel('I/bit');
% legend('0dB-IT','0dB-NP');
%% 恒模 正确
clc;
clear;
close all;

d =1;%信噪比
N = 128;%时间带宽积
p1=0:0.00001:1;
p0=1-p1;
n = 0 : 1 : N-1;
deltax = 1;
x = 0 : deltax : N-deltax;
% IP PFA PD
for kx = 1 : length(x)
    E(kx)=sqrt(1./(2*d*(N-x(kx))+N)).*exp(-(d*(N-x(kx)))^2./(2*(N+2*d*(N-x(kx)))));
    B(kx)=sqrt(1/(N+2*d*(N-x(kx))));
    C(kx)=exp(-(d*(N-x(kx)))^2./(2*N));
end
y=sum(E);
v=sum(B);
u=sum(C);
  pfa=(p1.*y)./(sqrt(N).*p0+p1.*y);
% PD=(p1.*v*sqrt(N))./(p1.*v*sqrt(N)+(1-p1).*u);
%NP PFA PD
z=zeros(1,100001);
for kx = 1 : length(x)
    z=z+1/2*erfc((sqrt(2*N)*erfcinv(2*pfa)-d*(N-x(kx)))/(sqrt(2*(2*(N-x(kx))*d+N))));
end
PD1=z./N;
HV=-p0.*log(p0)-p1.*log(p1);

p1_z0=pfa;%PFA
%p1_z0 = (p1.*y)./(sqrt(N).*p0+p1.*y);
p0_z0 = 1 - p1_z0;
hv_z0=-(p1_z0 .* log(p1_z0))-(p0_z0 .* log(p0_z0));

% p1_z1=PD;%PD
p1_z1 = (p1.*v*sqrt(N))./(p1.*v*sqrt(N)+p0*u);
p0_z1 = 1 - p1_z1;  
hv_z1 = -(p1_z1 .* log(p1_z1))-(p0_z1 .* log(p0_z1)) ;

IIT=HV-hv_z0-hv_z1;
semilogx(p1,IIT,'g');%IT
set(gca,'XScale','log')
hold on;

%NP
A=p1.*PD1./(p1.*PD1-p1.*pfa+pfa);
D=(1-p1).*(1-pfa)./(1-p1.*PD1+p1.*pfa-pfa);
INP=HV-(p1.*PD1-p1.*pfa+pfa).*(-A.*log(A)-(1-A).*log(1-A))-(1-p1.*PD1+p1.*pfa-pfa).*(-D.*log(D)-(1-D).*log(1-D));
semilogx(p1,INP,'r');%np
set(gca,'XScale','log')
xlabel('PFA'); ylabel('I/bit');
legend('0dB-IT','0dB-NP');

% %%%llkk
% N = 128;
% delta = 0.0001;
% P1= 0 : delta :1;
% SNR =0;
% x = 0:N-1;
% f = 0;
% g = 0;
% r = 0;
% NPP_D = 0;
% rho_2 = 10^(SNR/10);
% 
% % IT检测信息
% for i = 1:length(x)
%     f = f+sqrt(1/(2*rho_2*(N-x(i))+N)).*exp(-(rho_2*(N-x(i)))^2 ...
%     ./(4*rho_2*(N-x(i))+2*N));% 虚警概率检测统计量的分子
%     g = g+sqrt(1/(N+2*rho_2*(N-x(i))));%检测概率检测统计量的分子
%     r = r+sqrt(1/N)*exp(-rho_2^2*(N-x(i))^2/2/N);%检测概率检测统计量分母
% end
% temp1 = f/sqrt(N);
% ITP_FA = temp1*P1./(P1*temp1+1-P1);   
% temp2 = g/r;
% % temp2 = sqrt(2/pi)*(sqrt(N+2*N*rho_2)-sqrt(N+2*rho_2))/(erf(rho_2*sqrt(N/2))-erf(rho_2/sqrt(2*N)));%检测概率的统计量
% ITP_D = P1*temp2./(P1*temp2+1-P1);
% I = ones(1,10001);
% EP1Z0 = ITP_FA;
% EP1Z1 = ITP_D;
% EP0Z0 = I-ITP_FA;
% EP0Z1 = I-ITP_D;  
% % A = P1.*ITP_D./(P1.*ITP_D-P1.*P_FA+P_FA);
% % D = (1-P1).*(1-P_FA)./(1-P1.*ITP_D+P1.*P_FA-P_FA);
% % ITH_VV = (P1.*ITP_D-P1.*P_FA+P_FA).*(-A.*log2(A)-(1-A).*log2(1-A))+(1-P1.*ITP_D+P1.*P_FA-P_FA).*(-D.*log2(D)-(1-D).*log2(1-D));
% % IT_I = -P1.*log2(P1)-(1-P1).*log2(1-P1)-ITH_VV;
% H_V_Z = -EP1Z0.*log2(EP1Z0)-EP1Z1.*log2(EP1Z1)-EP0Z0.*log2(EP0Z0)-EP0Z1.*log2(EP0Z1);
% IT_I =-P1.*log2(P1)-(1-P1).*log2(1-P1)-H_V_Z;
% semilogx(P1,IT_I,'-','LineWidth',1)
% hold on;
% % NP检测信息
% P_FA = P1;
% for i = 1:length(x)
%     NPP_D = NPP_D+0.5*erfc((sqrt(2*N)*erfcinv(2*P_FA)-rho_2*(N-x(i)))/(sqrt(2*N+4*(N-x(i))*rho_2)))/N;
% end
% A = P1.*NPP_D./(P1.*NPP_D-P1.*P_FA+P_FA);
% D = (1-P1).*(1-P_FA)./(1-P1.*NPP_D+P1.*P_FA-P_FA);
% H_VV = (P1.*NPP_D-P1.*P_FA+P_FA).*(-A.*log2(A)-(1-A).*log2(1-A))...
%     +(1-P1.*NPP_D+P1.*P_FA-P_FA).*(-D.*log2(D)-(1-D  ).*log2(1-D));
% NP_I = -P1.*log2(P1)-(1-P1).*log2(1-P1)-H_VV;
% semilogx(P1,NP_I,'--','LineWidth',1)
% 
% % % z=zeros(1,10001);
% % % PFA=P_FA;
% % % d = rho_2;
% % % p1 = P1;
% % % for kx = 1 : length(x)
% % %     z=z+1/2*erfc((sqrt(2*N)*erfcinv(2*PFA)-d*(N-x(kx)))/(sqrt(2*(2*(N-x(kx))*d+N))));
% % % end
% % % PD1=z./N;A=p1.*PD1./(p1.*PD1-p1.*PFA+PFA);
% % % D=(1-p1).*(1-PFA)./(1-p1.*PD1+p1.*PFA-PFA);
% % % INP=-P1.*log2(P1)-(1-P1).*log2(1-P1)-(p1.*PD1-p1.*PFA+PFA).*(-A.*log(A)-(1-A).*log(1-A))-(1-p1.*PD1+p1.*PFA-PFA).*(-D.*log(D)-(1-D).*log(1-D));
% % % semilogx(PFA,INP,'r');%np
% 
% xlabel('先验概率P(1)');
% ylabel('检测信息I');                                                                                                                   
% legend("0dB-IT","0dB-NP")
% % axis([1e-3 1 -0.8 1])




%% 复高斯
% d =1;%信噪比
% N = 1024;%时间带宽积
% p1=0:0.0001:1;
% p0=1-p1;
% n = 0 : 1 : N-1;
% deltax = 1;
% x = 0 : deltax : N-deltax;
% % IP PFA PD
% for kx = 1 : length(x)
%     y(kx)=sqrt(1./((N-x(kx)*d^2)+2*d*(N-x(kx))+N)).*exp(-(d*(N-x(kx)))^2./(2*((N-x(kx)*d^2)+2*d*(N-x(kx))+N)));
%     v(kx)=1/sqrt((N-x(kx)*d^2)+2*d*(N-x(kx))+N);
%     u(kx)=exp((-d*(N-x(kx)))^2./(2*N));
% end
% y=sum(y);
% v=sum(v);
% u=sum(u);
% pfa=(p1.*y)./(sqrt(N).*p0+p1.*y);
% % PD=(p1.*v*sqrt(N))./(p1.*v*sqrt(N)+(1-p1).*u);
% %NP PFA PD
% z=zeros(1,10001);
% for kx = 1 : length(x)
%     z=z+1/2*erfc((sqrt(2*N)*erfcinv(2*pfa)-d*(N-x(kx)))./(sqrt(2*(2*(N-x(kx))*d^2+N+(N-x(kx))*d))));
% end
% PD1=z./N;
% HV=-p0.*log(p0)-p1.*log(p1);
% 
% % p1_z0=pfa;%PFA
% p1_z0 = (p1.*y)./(sqrt(N).*p0+p1.*y);
% p0_z0 = 1 - p1_z0;
% hv_z0=-(p1_z0 .* log(p1_z0))-(p0_z0 .* log(p0_z0));
% 
% % p1_z1=PD;%PD
% p1_z1 = (p1.*v*sqrt(N)/u)./(p1.*v*sqrt(N)/u+p0);
% p0_z1 = 1 - p1_z1;  
% hv_z1 = -(p1_z1 .* log(p1_z1))-(p0_z1 .* log(p0_z1)) ;
% 
% IIT=HV-hv_z0-hv_z1;
% semilogx(p1,IIT,'g');%IT
% set(gca,'XScale','log')
% hold on;
% 
% %NP
% A=p1.*PD1./(p1.*PD1-p1.*pfa+pfa);
% D=(1-p1).*(1-pfa)./(1-p1.*PD1+p1.*pfa-pfa);
% INP=HV-(p1.*PD1-p1.*pfa+pfa).*(-A.*log(A)-(1-A).*log(1-A))-(1-p1.*PD1+p1.*pfa-pfa).*(-D.*log(D)-(1-D).*log(1-D));
% semilogx(p1,INP,'r');%np
% set(gca,'XScale','log')
% xlabel('PFA'); ylabel('I/bit');
% legend('0dB-IT','0dB-NP');

