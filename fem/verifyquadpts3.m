function error = verifyquadpts3(n)
%% VERIFYQUADPTS3 examples and verfication on quadrature rules in 3-D.
%
%  error = verifyquadpts3(n) computes the error of n-th order quadrature
%  rule in a triangle. This is an example on the usage of quadrature points
%  and verification of qudarture order for approximating integrals in a
%  triangle.
%
% See also quadpts3, verifyquadpts
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if (n>5), n=5; end
% a reference tetrahedron
node = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
elem = [1 2 3 4];
volume = 1/6;
% get quadrature points
[lambda,weight] = quadpts3(n);
nQuad = size(lambda,1);
t1 = 0; 
t2 = 0;
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
         + lambda(p,2)*node(elem(:,2),:) ...
         + lambda(p,3)*node(elem(:,3),:) ...
         + lambda(p,4)*node(elem(:,4),:);
    t1 = t1 + weight(p)*f1(pxyz(1),pxyz(2),pxyz(3),n);
    t2 = t2 + weight(p)*f2(pxyz(1),pxyz(2),pxyz(3));
end                
t1 = t1*volume;
t2 = t2*volume;
error(1) = abs(t1 - 3/((n+1)*(n+2)*(n+3)));
error(2) = abs(t2 - (sin(1) + 1/2*cos(1)-1));
end

function z = f1(x,y,z,n)
z = x.^n + y.^n + z.^n;
end

function z = f2(x,y,z)
z = sin(x+y+z);
end

%% Results
% Let T be the triangle formed by (0,0), (1,0), and (0,1). 
%
% Error1 is for the integral 
% $\int _{T} x^n + y^n \, dxdy$. 
% It should be numerically exact. 
%
% Error2 is for the integral 
% $\int _{T} \sin(x+y) \, dxdy$.
% It decays as n increas.
%
% See the doc for qudrature rules in <matlab:ifemdoc('quadpts') quadpts>.
