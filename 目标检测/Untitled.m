a=0.00969;b=zeros(257,1);
size=length(b);
for i=size:-1:1
    a=a-0.000035;
    b(i,1)=a;
end