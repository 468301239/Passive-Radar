%%%ITµÄ¼ì²âÐÅÏ¢
clear all;close all;clc;
N=128;snr=[1 3.16];
x=0:0.000001:1;%%p(1)

for i=1:length(snr)
    d=snr(i);
    H_V=-x.*log2(x)-(1-x).*log2(1-x);
    
    Ez0_p0z0=(1-x)./(1-x+x./sqrt(1+2*d).*exp((-N*d^2)/(2*(1+2*d))));
    Ez0_p1z0=(x./sqrt(1+2*d).*exp((-N*d^2)/(2+4*d)))./(1-x+x./sqrt(1+2*d).*exp((-N*d^2)/(2*(1+2*d))));
    Ez1_p0z1=(1-x)./(1-x+x./sqrt(1+2*d)*exp(0.5*N*d^2));
    Ez1_p1z1=x./(x+sqrt(1+2*d).*(1-x).*exp(-0.5*N*d^2));
    
    H_VZ=-Ez0_p0z0.*log2(Ez0_p0z0)-Ez0_p1z0.*log2(Ez0_p1z0)-Ez1_p0z1.*log2(Ez1_p0z1)-Ez1_p1z1.*log2(Ez1_p1z1);
    I=H_V-H_VZ;
    semilogx(x,I);
    hold on;
end
legend('SNR=0dB','SNR=5dB');