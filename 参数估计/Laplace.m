%% 后验的分子绘图
clc;clear;close all;
k=1.380649*10^-23;T=273.15+25;Pi=3.1415926;sigma=1;N0=1;%N0=k*T;
N=16;x0=N/2;n=0:N;SNR=10;
alpha=sqrt(SNR*N0);
% wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% w=sum(wn);

mu=0;                      %均值
sigma=1;                 %标准差，方差的开平方
b=sigma/sqrt(2);      %根据标准差求相应的b
aa=rand(1,10000)-0.5;    %生成(-0.5,0.5)区间内均匀分布的随机数列;
wn=mu-b*sign(aa).*log(1-2*abs(aa)); %生成符合拉普拉斯分布的随机数列
w=sum(wn);
Mean=mean(wn);
Std=std(wn);
% hist(wn,100);

x=0:0.1:N-0.1;
phi=Pi/2;phi0=Pi/2;
in=0;
for n=1:N
    qq=Mean/alpha;
    in=in+abs(exp(1i*phi)*sinc(n-x0)-exp(1i*phi).*sinc(n-x)+qq);
end
result=besseli(0,alpha/sigma.*in);
plot(x,result,'r','LineWidth',2);

% %% 贝塞尔函数里sinc-sinc
% clc;clear;close all;
% N=100;
% x0=N/4;
% x=0:0.1:N;
% y=0;
% for n=0:N
%     y=y+besseli(0,sinc(n-x0)-sinc(n-x));
% end
% in=@(x1) besseli(0,sinc(n-x0)-sinc(n-x1));
% A=integral(in,0,N);
% plot(x,y/A);

% %% sinc-sinc
% clc;clear;close all;
% N=20;
% x0=N/2;
% x=0:0.1:N;
% y=0;
% for n=0:N
%     y=y+sinc(n-x0)-sinc(n-x);
% end
% plot(x,y);

% %% 后验的分子绘图 x.x0合并后
% clc;clear;close all;
% k=1.380649*10^-23;T=273.15+25;Pi=3.1415926;N0=1;%N0=k*T;
% N=8;x0=N/2;n=0:N;SNR=10;
% alpha=sqrt(SNR*N0);
% % wn=sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
% % w=sum(wn);

% mu=0;                      %均值
% sigma=1;                 %标准差，方差的开平方
% b=sigma/sqrt(2);      %根据标准差求相应的b
% aa=rand(1,10000)-0.5;    %生成(-0.5,0.5)区间内均匀分布的随机数列;
% wn=mu-b*sign(aa).*log(1-2*abs(aa)); %生成符合拉普拉斯分布的随机数列
% w=mean(wn);
% 
% 
% x=0:0.1:N;
% phi=Pi/2;phi0=Pi/2;
% in=0;
% for n=1:N
%     display(n);
%     in=in+abs(exp(1i*phi)*(sin(n-x0)/(n-x0)^2-cos(n-x0)/(n-x0)).*abs(x-x0)+w/alpha);
% end
% result=besseli(0,alpha/sigma.*in);
% plot(x,result,'r','LineWidth',2);

