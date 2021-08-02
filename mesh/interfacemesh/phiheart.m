function z = phiheart(p)
%% PHIHEART level set function for a heart
%
% Author: Huayi Wei <weihuayi@xtu.edu.cn>

x = p(:,1) - 0.02; y = p(:,2) -0.02;
z = (x.^2 + y.^2 - 0.38).^3 - x.^2.*y.^3;
end