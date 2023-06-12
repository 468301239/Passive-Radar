clear all;close all;clc;

%% 卡方加高斯多图

sw=1;
a=2;
M=[2,2^2,2^3,2^4];
pi=3.1416;
z=-10:1:100;
%P=sigma^(T-2)*(alpha*T)*(T/2-1)/(4*gamma(T)*sqrt(pi)).*exp(-z.^2./(4*alpha*T*sigma^2)).*(sigma*sqrt(alpha*T)*gamma(T/2).*hypergeom(T/2,1/2,(z-alpha*T*sigma^2).^2./4*alpha*T*sigma^2)+(z-alpha*T*sigma^2).*gamma((T+1)/2)).*hypergeom((T+1)/2,3/2,(z-alpha*T*sigma^2).^2./4*alpha*T*sigma^2);
figure;
for ii=1:length(M)
    N=M(ii);
    P=(a*N)*(N/2-1)/(4*gamma(N)*sqrt(pi)).*exp(-z.^2./(4*a*N)).*(sqrt(a*N)*gamma(N/2).*hypergeom(N/2,1/2,(z-a*N).^2./(4*a*N))+(z-a*N).*gamma((N+1)/2).*hypergeom((N+1)/2,3/2,(z-a*N).^2./(4*a*N*sw^2)));
%     P=(sw^-(2*N+1)*a^(N/2-3)*exp(-z.^2./(8*a^2*N*sw^2))*(1/(sw^2*a^2*N)).^(-N/2))./(gamma(N)*a*sqrt(pi*N)).*(sqrt(2)*gamma(N/2)*hypergeom(N/2,1/2,(z-a*N*a^2).^2./(8*a^2*N*sw^2))+(z-a*N*a^2).*gamma((N+1)/2).*sqrt(1/(N*sw^2*a^2)).*hypergeom((N+1)/2,3/2,(z-a*N*a^2).^2./(8*a^2*N*sw^2)));
    %此为式1.11
    hold on;plot(z,P,'LineWidth',1);grid on;
% %    Q=1/(sqrt(2*pi*(2*a*N+2*2*N)))*exp(-(z-(2*N+1)).^2./(2*(2*a*N+2*2*N)));
% %    Q=1/(sqrt(2*pi*(2*a*N*sw^2+2*2*N)))*exp(-(z-(2*N-1)).^2./(2*(4*a^2*N*sw^2+2*2*N)));
       Q=1/(sqrt(2*pi*(4*a^2*N*sw^2+2*2*N)))*exp(-(z-(2*N-1+(a^2.*N))).^2./(2*(4*a^2*N*sw^2+2*2*N)));
    hold on;plot(z,Q,'r');grid on;
end
%legend('2','4','8','16','100')
xlabel('x');
ylabel('PDF');
legend('N=2','N=4','N=8','N=16')

%% N=10时卡方加高斯


%% 验证单独的卡方和高斯_正确
% N=4;
% sw=1;% sw=1;%sigma_w
% a=2;%alpha
% pi=3.1416;
% z=-10:1:100;
% %卡方自由度k=2N
% N=12;%窗长
% %P=(a*N)^(N/2-1)/(4*gamma(N)*sqrt(pi)).*exp(-z.^2./(4*a*N)).*(sqrt(a*N)*gamma(N/2).*hypergeom(N/2,1/2,(z-a*N).^2./(4*a*N))+(z-a*N).*gamma((N+1)/2).*hypergeom((N+1)/2,3/2,(z-a*N).^2./(4*a*N)));
% %P=sw^(N-2)*(a*N)^(N/2-1)/(4*gamma(N)*sqrt(pi)).*exp(-z.^2./(4*a*N*sw^2)).*(sw*sqrt(a*N)*gamma(N/2).*hypergeom(N/2,1/2,(z-a*N*sw^2).^2./(4*a*N*sw^2))+(z-a*N*sw^2).*gamma((N+1)/2).*hypergeom((N+1)/2,3/2,(z-a*N*sw^2).^2./(4*a*N*sw^2)));
% % P=2^(N/2-3)*(a*sw)^(N-2)*N^(N/2-1).*exp(-z.^2./(8*a^2*N*sw^2)).*(a*sw*sqrt(2*N)*gamma(N/2).*hypergeom(N/2,1/2,(z-2*a^2*N*sw^2).^2./(8*a^2*N*sw^2))+(z-2*a^2*N*sw^2).*gamma((N+1)/2).*hypergeom((N+1)/2,3/2,(z-2*a^2*N*sw^2).^2./(8*a^2*N*sw^2)));
% P=(sw^(-2*N-1)*a^(N/2-3)*exp(-z.^2./(8*a^2*N*sw^2))*(1/(sw^2*a^2*N)).^(-N/2))./(gamma(N)*a*sqrt(pi*N)).*(sqrt(2)*gamma(N/2)*hypergeom(N/2,1/2,(z-a*N*a^2).^2./(8*a^2*N*sw^2))+(z-a*N*a^2).*gamma((N+1)/2).*sqrt(1/(N*sw^2*a^2)).*hypergeom((N+1)/2,3/2,(z-a*N*a^2).^2./(8*a^2*N*sw^2)));
% figure;
% %hold on;plot(z,P);grid on;
% 
% %Q=1/(sqrt(2*pi*(2*a*N+2*2*N)))*exp(-(z-(2*N-1)).^2./(2*(2*a*N+2*2*N)));
% %Q=1/(sqrt(2*pi*(2*a*N*sw^2+2*2*N)))*exp(-(z-(2*N-1)).^2./(2*(2*a*N*sw^2+2*2*N)));
% Q=1/(sqrt(2*pi*(4*a^2*N*sw^2+2*2*N)))*exp(-(z-(2*N-1+(a^2.*N))).^2./(2*(4*a^2*N*sw^2+2*2*N)));
% 
% hold on;plot(z,Q);grid on;
% xlabel('x');
% ylabel('PDF（N）');
% legend('chisquare+gauss','gauss');

% a=2;
% x=0:100;
% figure;
% P=x.^(N-1).*exp(-x./2)/(2^N*gamma(N));
% hold on;plot(x,P);grid on;
% Q=1/(2*sw*sqrt(pi*a*N)).*exp(-x.^2./(4*a*N*sw^2));
% hold on;plot(x,Q);grid on;




%% 1F1对应hypergeom没问题
% T=4;
% a=2;
% x=0.1:0.1:10;
% Q=hypergeom(2,1,x);
% P=x.^(-2).*hypergeom([2,2],[],-1./x);
% P=gamma(0)/gamma(2).*hypergeom(2,1,x)+gamma(0)/gamma(2).*x.^0.*hypergeom(2,1,x);
% plot(x,P,'r','LineWidth',1)
