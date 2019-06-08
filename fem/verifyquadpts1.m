function error = verifyquadpts1(n)
%% VERIFYQUADPTS1 examples and verfication on quadrature rules in 1-D.
%
%  error = verifyquadpts1(n) computes the error of n-th order quadrature
%  rule in a triangle. This is an example on the usage of quadrature points
%  and verification of qudarture order for approximating integrals in a
%  triangle.
%
% See also quadpts
%
% Added by Huayi Wei
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if (n>10), n=10; end
% a reference triangle

a = 0;
b = 1;
l = 1;


% get quadrature points
[lambda,weight] = quadpts1(n);
nQuad = size(lambda,1);
t1 = 0; 
t2 = 0;
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    px = lambda(p,1)*a+ lambda(p,2)*b;

    t1 = t1 + weight(p)*f1(px,n);
    t2 = t2 + weight(p)*f2(px);
end                
t1 = t1*l;
t2 = t2*l;
error(1) = abs(t1 - 1/(2*n));
error(2) = abs(t2 - (- cos(1)-1));
end

function z = f1(x, n)
z = x.^(2*n-1);
end

function z = f2(x)
z = sin(x);
end

