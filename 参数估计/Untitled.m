clc;clear;close all; 
N0=1;
N=20000;
wn=sqrt(N0/2)*(randn(1,N)+1i*randn(1,N));
 m1=abs(mean(wn))
 s1=std(wn)