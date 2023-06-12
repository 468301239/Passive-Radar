%%%IT  NPµÄ¼ì²âÐÅÏ¢
clear all;close all;clc;
N=128;snr=[5];
x=0:0.0001:1;%%p(1)

%%%IT
for i=1:length(snr)
    d=10^(snr(i)/10);
    H_V=-x.*log2(x)-(1-x).*log2(1-x);
    
    Ez0_p0z0=(1-x)./(1-x+x./sqrt(1+2*d).*exp((-N*d^2)/(2*(1+2*d))));
    Ez0_p1z0=(x./sqrt(1+2*d).*exp((-N*d^2)/(2+4*d)))./(1-x+x./sqrt(1+2*d).*exp((-N*d^2)/(2*(1+2*d))));
    Ez1_p0z1=(1-x)./(1-x+x./sqrt(1+2*d)*exp(0.5*N*d^2));
    Ez1_p1z1=x./(x+sqrt(1+2*d).*(1-x).*exp(-0.5*N*d^2));
    
    H_VZ=-Ez0_p0z0.*log2(Ez0_p0z0)-Ez0_p1z0.*log2(Ez0_p1z0)-Ez1_p0z1.*log2(Ez1_p0z1)-Ez1_p1z1.*log2(Ez1_p1z1);
    I=H_V-H_VZ;
    semilogx(x,I,'r','LineWidth',1);ylim([0 1]);
    hold on;
end

%%%NP
for i=1:length(snr)
    d=10^(snr(i)/10);
    H_V=-x.*log2(x)-(1-x).*log2(1-x);
    
    PFA=(x./sqrt(1+2*d)*exp((-N*d^2)/(2+4*d)))./(1-x+x./sqrt(1+2*d)*exp((-N*d^2)/(2+4*d)));
    PD=0.5*erfc((sqrt(2*N)*erfcinv(2*PFA)-N*d)./(sqrt(2*N+4*N*d^2)));
    A=(x.*PD)./(x.*PD-x.*PFA+PFA);
    D=((1-x).*(1-PFA))./(1-x.*PD+x.*PFA-PFA);
    I=H_V-(x.*PD-x.*PFA+PFA).*(-A.*log2(A)-(1-A).*log2(1-A))-(1-x.*PD+x.*PFA-PFA).*(-D.*log2(D)-(1-D).*log2(1-D));
    semilogx(x,I,'b--','LineWidth',1);ylim([0 1]);
    hold on;
end
xlabel('P_1');
ylabel('I');
% legend('IT','NP');