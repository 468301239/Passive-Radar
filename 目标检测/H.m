% function [hx] = Hx(a)
% % hx=0;
%     Length=length(a);
%     for i=1:1:Length
%         if i==1
%             b=0;
%         end
%         b=b-a(i)*log2(a(i));
%     end
%     hx=b;
% end

%POWERED BY ZHUZHENGTAO
function H = H(a)
    I=log2(a);
    H=-I*a';
end

% function [hx] = Hx(a)
% %求信息熵
% %   此处显示详细说明
%     b=arrayfun(@plnp,a);
%     hx=sum(b);
% end
% function [y] = plnp(p)
% %y=-p*log2(p)
%     y=-p*log2(p);
% end



