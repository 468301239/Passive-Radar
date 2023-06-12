N=4;k=2*N;sigma_w=1;alpha=2;
syms x;
y=pw_x(x).*psw_x(z-x);
pr_z=int(y,0,inf);
figure;
for z=-20:0.1:60
    pr_z;
    hold on;plot(z,pr_z);grid on;
end

function y=pw_x(x)
    N=4;k=2*N;sigma_w=1;alpha=2;
    y=(x.^(k-1)).*exp(-x./(2.*sigma_w.^2))...
        ./(2.^(k/2).*sigma_w.^k.*gamma(k/2));
end
function y=psw_x(x)
    N=4;k=2*N;sigma_w=1;alpha=2;
    y=exp(-x.^2./(2*4*alpha.^2*N*sigma_w.^2));
end