function error = verifyquadquadpts(n)
%% VERIFYQUADQUADPTS verify the quadrature formula on [0,1]^2
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (n>10), n=10; end

a = 0;
b = 1;

c = 0;
d = 1;

l = 1;


% get quadrature points
[pts,weight] = quadquadpts(n);
nQuad = size(pts,1);
t1 = 0; 
t2 = 0;
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    px = pts(p,1)*a + (1-pts(p,1))*b;
    py = pts(p,2)*c + (1-pts(p,2))*d;

    t1 = t1 + weight(p)*f1(px,py, n);
    t2 = t2 + weight(p)*f2(px,py);
end                
t1 = t1*l;
t2 = t2*l;
error(1) = abs(t1 - 1/(4*n*n));
error(2) = abs(t2 - (- cos(1)+1)*(-cos(1)+1));
end

function z = f1(x, y, n)
z = x.^(2*n-1).*y.^(2*n-1);
end

function z = f2(x,y)
z = sin(x)*sin(y);
end