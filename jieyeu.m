clc;clear;close all;

N=2^4;
n=-N/2:1:N/2-1;
x=-N/2:1:N/2-1;
x0=0;
% y=sum(stepfun(n,x).*stepfun(n,x0));
for i=1:length(x)
    y(i)=sum(real(n>x(i)).*real(n>x0));
end
figure;
plot(x,y,'LineWidth',2);