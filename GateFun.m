N=2^4;
n=0:1:N-1;
x=0:1:N-1;
x0=2;
g=heaviside(n+0.4)-heaviside(n-0.4);
% for i=1:length(x)
% y(i)=sum(g.*(heaviside(n+0.4-x(i))-heaviside(n-0.4-x(i))));
% end
% plot(x,y);
yy=heaviside(n+0.01-x0)-heaviside(n-0.01-x0);
yyy=(heaviside(n+0.01-x0)-heaviside(n-0.01-x0)).^2;
y=sum((heaviside(n+0.01-x0)-heaviside(n-0.01-x0)).^2);

%%
N=2^5;
n=-N/2:1:N/2-1;
x=-N/2:0.1:N/2-1;
x0=0;
g=heaviside(n+0.4)-heaviside(n-0.4);
for i=1:length(x)
% y(i)=sum(g.*(heaviside(n+0.4-x(i))-heaviside(n-0.4-x(i))));
y(i)=sum(men(N,x0).*men(N,x(i)));
end
plot(x,y);

