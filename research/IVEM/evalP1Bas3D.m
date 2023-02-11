function f = evalP1Bas3D(bas, x, y, z, dind)

%% USAGE: evaluate P1 solution or interpolation at certain point(s).
%
% INPUTS:
% bas --- coefficient of basis [c1,c2,c3,c4], ci can be vectors
%            c1 + c2*x + c3*y + c4*z
% x, y, z --- coordinates of the points to evaluate can becolumn vectors,
% dind --- derivative info for FE function
%            d = [0,0,0]: function value
%            d = [1,0,0]: Dx value
%            d = [0,1,0]: Dy value
%            d = [0,0,1]: Dz value
% OUTPUTS:
% f --- the value of the FE/DG solution or interpolation at point(x,y). 
%
% Last Modified: 06/19/2020 by Xu Zhang
% Last Modified: 07/02/2020 by Xu Zhang
%%
one = ones(size(x));
c1 = bas(:,1); c2 = bas(:,2); c3 = bas(:,3); c4 = bas(:,4);

if dind(1) == 0 && dind(2) == 0 && dind(3) == 0
    f = c1.*one + c2.*x + c3.*y + c4.*z;
elseif dind(1) == 1 && dind(2) == 0 && dind(3) == 0
    f = c2.*one;
elseif dind(1) == 0 && dind(2) == 1 && dind(3) == 0
    f = c3.*one;
elseif dind(1) == 0 && dind(2) == 0 && dind(3) == 1
    f = c4.*one;
elseif sum(dind) > 1
    f = zeros(size(x));
end