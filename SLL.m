clc;
clear;
close all;

num=50;
N=64;
n=-N/2:N/2-1;
t=2;
x=-N/2:N/2-t;
SNR=[-30,-20,-10,0];
p1=0:1e-4:1;

temp1_z0 = zeros(1,N);
pfa=p1;
plot(p1,pfa,'k-','LineWidth',1);hold on;
for snr =1:length(SNR)
    alpha=1;
    rho2 = 10^(SNR(snr)/10);
    N0=alpha^2/rho2;
    for run = 1:num
        w=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));        
        for kx = 1:length(x)
            
            temp1_z0(kx) = exp(-rho2)*besseli(0,2*alpha/N0*abs(mean((heaviside(n-x(kx))-...
                heaviside(n-x(kx)-t)).*w)));
        end
        gamma1_z0(run)=1/N*sum(temp1_z0);
    end
    P_FA = p1*mean(gamma1_z0)./(1-p1+p1*mean(gamma1_z0));
    plot(p1,P_FA,'--','LineWidth',1);hold on;
end

legend('P_F_A=P(1)','SNR=-30','SNR=-20','SNR=-10','SNR=0')
xlabel('P(1)');
ylabel('P_F_A');