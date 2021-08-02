function [z,n] = phiflower(p)
%% PHIFLOWER level set function of a flower
%
% Author: Huayi Wei <weihuayi@xtu.edu.cn>

petal = 5;
xc = 0;
yc = 0;
r0 = 0.5;
r1 = 0.2;
theta = atan2((p(:,2) - yc),(p(:,1) - xc));
r = r0 + r1*sin(petal*theta);

z = (p(:,2)-yc).^2 + (p(:,1)-xc).^2 - r.^2;

if nargout == 2
    rp = petal*r1*cos(petal*theta);
    n1 = rp.*cos(theta) - r.*sin(theta);
    n2 = rp.*sin(theta) + r.*cos(theta);
    l = sqrt(n1.^2+n2.^2);
    n = [n2./l, -n1./l];
end
end