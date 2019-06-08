function [u,p,err] = uzawa(A,B,f,g,r,omega)
%% UZAWA augmented Uzawa method for solving Stokes equations
% 
%      |A  B'| |u|  = |f|
%      |B  0 | |p|  = |0|
%
% The parameter r is used in the augmented uzawa method. 
%
%      |A+rB'B  B'| |u|  = |f|
%      |B       0 | |p|  = |0|
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin<=4
    r = 0;    omega = 2;
elseif nargin <=5;
    omega = max(r+1,2);
end

Np = size(B,1);
p = zeros(Np,1);
u = zeros(length(f),1);

A = A + r*(B'*B);
f = f + r*(B'*g);
err = 1;
k = 1;
while (err>1e-6) && (k<10000)
    uold = u;
    u = A\(f-B'*p);
    p = p + omega*(B*u-g);
    err = sqrt((u-uold)'*A*(u-uold))/sqrt(u'*A*u);
    fprintf('#dof: %8.0u, Uzawa iter: %2.0u, err = %12.8g\n \n',...
         length(u)+length(p), k, err);
    k = k+1;
end