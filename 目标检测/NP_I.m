%%%NPµÄI
clear all;close all;clc;
N=128;snr=[3.16];
x=0:0.0001:1;%%p(1)

for i=1:length(snr)
    d=snr(i);
    H_V=-x.*log2(x)-(1-x).*log2(1-x);
    
    PFA=(x./sqrt(1+2*d)*exp((-N*d^2)/(2+4*d)))./(1-x+x./sqrt(1+2*d)*exp((-N*d^2)/(2+4*d)));
    PD=0.5*erfc((sqrt(2*N)*erfcinv(2*PFA)-N*d)./(sqrt(2*N+4*N*d^2)));
    A=(x.*PD)./(x.*PD-x.*PFA+PFA);
    D=((1-x).*(1-PFA))./(1-x.*PD+x.*PFA-PFA);
    I=H_V-(x.*PD-x.*PFA+PFA).*(-A.*log2(A)-(1-A).*log2(1-A))-(1-x.*PD+x.*PFA-PFA).*(-D.*log2(D)-(1-D).*log2(1-D));
    
    semilogx(x,I,'r');
    hold on;
end
%legend('SNR=0dB','SNR=5dB');