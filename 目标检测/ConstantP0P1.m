% �࿨����˹�������һ��
% clear all;close all;clc;
% sw=1;%Ϊʲô�����Ӱ���˹���ƵĽ����
% % sw=1ʱ�����ֲ��ļ�ֵ����k-2,sw=!1ʱ��ֵ����(k-2)*sw^2
% x=0:50;
% 
% figure;
% 
% k=3:20;
% for ii=1:length(k)
%     P=x.^(k(ii)/2-1).*exp(-x/(2*sw^2))./(2^(k(ii)/2)*sw^k(ii)*gamma(k(ii)/2));
%     hold on;plot(x,P);grid on;
%     
% end
% 
% k1=20;
% max=(k1-2)^(k1/2-1)*exp(-(k1-2)/(2*sw^2))/(2^(k1/2)*sw^k1*gamma(k1/2));
% sigma=1/(sqrt(2*pi)*max);
% P0=1/(sqrt(2*pi)*sigma)*exp(-(x-k1+2).^2/(2*sigma^2));
% hold on;plot(x,P0,'r','LineWidth',1);grid on;
% 
% xlabel('x');
% ylabel('ChiSquare');
% legend('k=20-30');




% �����ֲ���һ��ͼ��
clear all;close all;clc;

sw=1;
z=-20:0.1:60;

% figure;

% k=20:30;
% for ii=1:length(k)
%     P=x.^(k(ii)/2-1).*exp(-x/2)./(2^(k(ii)/2)*gamma(k(ii)/2));
%     hold on;plot(x,P);grid on;
%     
% end

N=4;
k1=2*N;
max=(k1-2)^(k1/2-1)*exp(-(k1-2)/(2*sw^2))/(2^(k1/2)*sw^k1*gamma(k1/2));

sigma=1/(sqrt(2*pi)*max);
P0=1/(sqrt(2*pi)*sigma).*exp(-(z-k1+2).^2./(2*sigma^2));
hold on;plot(z,P0,'LineWidth',1);grid on;


%ֱ���÷�����Ʋ��Զ
% sigma2=2*(k1-2);
% Q2=1/(sqrt(2*pi)*sigma2).*exp(-(x-k1+2).^2./(2*sigma2^2));
% hold on;plot(x,Q2,'LineWidth',1);grid on;

a=2;%alpha

mu=2*N-1+a^2*N;%���Ƶĸ�˹��ֵ
s=4*N+sw^2;%���Ƶĸ�˹����
P1=1/(sqrt(2*pi*(4*N+sw^2))).*exp(-(z-(2*N-1+a^2*N)).^2./(2*(4*N+sw^2)));
hold on;plot(z,P1,'LineWidth',1);grid on;

xlabel('x');
ylabel('ChiSquare');
legend('P0','P1');

solve(P0-P1,z)

% k=3:10ÿ��ͼһ���������˹�ıȽ�
% clear all;close all;clc;
% x=0:0.1:50;
% k=3:10;%���ɶ�ȡֵ��Χ
% pi=3.1416;
% for ii=1:length(k)
%     P=x.^(k(ii)/2-1).*exp(-x/2)./(2^(k(ii)/2)*gamma(k(ii)/2));
%     max=(k(ii)-2)^(k(ii)/2-1)*exp(-(k(ii)-2)/2)/(2^(k(ii)/2)*gamma(k(ii)/2));
%     sigma=1/(sqrt(2*pi)*max);
%     Q=1/(sqrt(2*pi)*sigma).*exp(-(x-k(ii)+2).^2./(2*sigma^2));
%     figure;
%     hold on;plot(x,P);grid on;
%     hold on;plot(x,Q,'r','LineWidth',1);grid on;
%     
% xlabel('x');
% ylabel('ChiSquare');
% end


