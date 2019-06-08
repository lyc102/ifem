function verifyquadpts
%% VERIFYQUADPTS examples and verfication on quadrature rules.
%
%  error = verifyquadpts(n) computes the error of n-th order quadrature
%  rule in a triangle. This is an example on the usage of quadrature points
%  and verification of qudarture order for approximating integrals in a
%  triangle.
%
% See also quadpts
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

% a reference triangle
node = [0 0; 1 0; 0 1];
elem = [1 2 3];
area = 0.5;
err = zeros(9,2);
for n = 1:9    
    % get quadrature points
    [lambda,weight] = quadpts(n);
    nQuad = size(lambda,1);
    t1 = 0; 
    t2 = 0;
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        t1 = t1 + weight(p)*f1(pxy(1),pxy(2),n);
        t2 = t2 + weight(p)*f2(pxy(1),pxy(2));
    end                
    t1 = t1*area;
    t2 = t2*area;
    err(n,1) = abs(t1 - 2/((n+1)*(n+2)));
    err(n,2) = abs(t2 - (sin(1) - cos(1)));
end
display('Table: Error of quadrature for two smooth functions')
colname = {'n','x^n + y^n','sin(x+y)'};
disptable(colname,(1:9)',[],err(:,1),'%0.5e',err(:,2),'%0.5e');
end

function z = f1(x,y,n)
z = x.^n + y.^n;
end

function z = f2(x,y)
z = sin(x+y);
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
