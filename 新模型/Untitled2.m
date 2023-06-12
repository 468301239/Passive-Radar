N0 = 1e-20;rho =0.01;alpha = sqrt(rho*N0);
P1 = 0:0.0001:1;P0 = 1-P1;
temp =1;
F = 0;
% num=100;
% temp = zeros(1,length(x))+1;%预置被积函数
% wn = zeros(N,num);
% F = zeros(num,length(x));%预置似然比
P_FA = P1; 
plot(P1,P_FA,"LineWidth",1)
hold on;
for N = [128 512 1024]
    n = 1:N;
    wn = sqrt(N0/2)*(randn(1,length(n))+1i*randn(1,length(n)));
    for x =1:N
        for kx = x: N
            temp=temp*(exp(-rho)*besseli(0,2*alpha/N0*abs(wn(kx)))); 
        end
        F = F+temp;     
    end
    P_FA2 = P1*F./(N*P0+P1*F);
    plot(P1,P_FA2,'--',"LineWidth",1)
    grid on;
    hold on;
    legend("P_F_A=P1","N=128","N=512","N=1024")
    xlabel('\pi(1)');
    ylabel('P_F_A');  
end