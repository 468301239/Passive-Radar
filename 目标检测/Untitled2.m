clc;clear;close all;
N=5e6;
deltax=0.01;                     
x=0:deltax:N;
T_1_z0=0;
SNR=1;
p1=0:0.01:1;
p0=1-p1;

for kx=1:length(x)
    T_1_z0=T_1_z0+1/N*sqrt(N/(N+2*(N-x(kx))*SNR))*exp(-(N-x(kx))^2*SNR^2./(2*N+4*(N-x(kx))*SNR))*deltax;
end
p0z0=p0./(p0+p1.*T_1_z0);