clear all;clc;

SNR = -15:2:25;
rho2 = 10.^(SNR/10);%信噪比
N = 2^4;
x10 = 0;
deltax = 0.05;
x1 = -N/2 : deltax : N/2-deltax;
n = -N/2 : 1 : N/2-1;
P1 = 1;
e=2.718;
R1 = sinc(n-x10)'*sinc(n-x10);

EEB=zeros(1,length(SNR));
Num =10;
MSE=zeros(1,length(SNR));
SIGMAEE=zeros(1,length(SNR));
SIGMAEE1=zeros(1,length(SNR));
SIGMAEE2=zeros(1,length(SNR));

tic
for snr=1:length(SNR)
    N0=2*P1^2/rho2(snr);
    CRB(snr)=1/(pi^2/3*rho2(snr));% 恒模的CRB
%   phi0=2*pi*rand(1,1);
    sapxh=zeros(1,Num);
    h=zeros(1,Num);
    row=zeros(1,Num);
    summ1=0;
    for run=1:Num        
        p_z_x = zeros(1,length(x1));        
        wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
        for kx=1:length(x1)
            p_z_x(1,kx)=2*pi*besseli(0,rho2(snr)*...
                        abs(sinc(x1(kx)-x10)+1/P1*...
                        sum(wn.*sinc(n-x1(kx)))));
        end     
        p_x_z = p_z_x./sum(p_z_x)/deltax; 
        h(1,run)=-sum(p_x_z.*log2(p_x_z+eps))*deltax; %理论熵，eps防止出现log0
        f=zeros(1,length(x1));
        summ=0;       
        for i = 1:length(x1)           
               summ=summ+p_x_z(i)*deltax;
               f(1,i)=summ;  %概率分布函数       
        end        
        randx = rand(1);
        for i = 1:length(x1)          
           if (randx<f(1,i))
             sapxh(1,run)=-log2(p_x_z(i));      %SAP熵                             
             break;
           end 
        end
        
        maxp=max(p_z_x);
        row(1,run)=find(p_z_x==maxp);
%         row(1,run)=find(p_z_x==maxp);
        Location=(row-1).*deltax-N/2;        
        summ1=summ1+Location(1,run)^2;            
    end
    x3=x1./deltax+N/2/deltax+1;
    p1=hist(row,x3);
    p4=p1./sum(p1)/deltax;    
    maxxh=-sum(p4.*log2(p4+eps))*deltax;  %MLE熵  
    SIGMAEE(1,snr)=2^(2*mean(h))/(2*pi*e);
    SIGMAEE1(1,snr)=2^(2*maxxh)/(2*pi*e);
    MSE(1,snr)=summ1/Num;
    SIGMAEE2(1,snr)=2^(2*mean(sapxh))/(2*pi*e);
end
toc

figure;semilogy(SNR,SIGMAEE,SNR,SIGMAEE1,SNR,MSE,SNR,CRB,SNR,SIGMAEE2)
xlabel('SNR');ylabel('ERROR');
legend('EE','MLE-EE','MLE-MSE','CRB','SAP-EE');