function [output] = gate(x0)

width=2;
% output=heaviside(x0+width/2)-heaviside(x0-width/2);

output=u(x0+width/2)-u(x0-width/2);
end

% function [output] = men(N,width,center)
% %MEN �˴���ʾ�йش˺�����ժҪ
% %   �˴���ʾ��ϸ˵��
% N=N;
% n=0:1:N-1;
% 
% output=heaviside(n+width/2-center)-heaviside(n-width/2-center);
% end