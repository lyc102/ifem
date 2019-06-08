function [phi, Dphi, lambda, J] = quadbasisNew(node,elem, GaussPt, elemType)
%% QUADBASIS the basis and their derivatives and the det of Jacobi.
% 
%   [phi, Dphi, J] = quadbasis(node,elem, GaussPts) compute the values of 
%   bilinear basis and their derivatives and Jacobi det at a gauss point GaussPt in a
%   quadrilateral.
%   elemType: 'Q1', 'Q2', ....
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

NT = size(elem,1);
xi = GaussPt(:,1);
eta = GaussPt(:,2);

%% basis function and their derivatives
%   [ eta - 1, xi - 1
%     1 - eta, -xi
%     eta,      xi
%     -eta,    1 - xi ]
% for a point (x, y) in real element,  there exist a point in reference
% element (xi, eta), satisfy
%    x = x_1*phi_1 + x_2*phi_2 + x_3*phi_3 + x_4 *phi_4
%    y = y_1*phi_1 + y_2*phi_2 + y_3*phi_3 + y_4 *phi_4
% So the Jacobi matrix is 
%
% [ x_2 - x_1 + eta*(x_1 - x_2 + x_3 - x_4), x_4 - x_1 + xi*(x_1 - x_2 + x_3 - x_4)
%   y_2 - y_1 + eta*(y_1 - y_2 + y_3 - y_4), y_4 - y_1 + xi*(y_1 - y_2 + y_3 - y_4)], 
%
% 
phi = zeros(4,1);
phi(1) = (1 - xi).*(1 - eta);
phi(2) = xi.*(1 - eta);
phi(3) = xi.*eta;
phi(4) = (1 - xi).*eta;

a = node(elem(:,1),1) - node(elem(:,2),1) + node(elem(:,3),1) - node(elem(:,4),1);
b = node(elem(:,1),2) - node(elem(:,2),2) + node(elem(:,3),2) - node(elem(:,4),2);

t1 = node(elem(:,2),1) - node(elem(:,1),1);
t2 = node(elem(:,4),1) - node(elem(:,1),1);

t3 = node(elem(:,2),2) - node(elem(:,1),2);
t4 = node(elem(:,4),2) - node(elem(:,1),2);


a11  = t1 + eta*a;
a12 = t2 + xi *a;
a21 = t3 + eta*b;
a22 = t4 + xi *b;

J = abs(a11.*a22 - a21.*a12);

Dphi = zeros(NT, 2, 4);

Dphi(1:NT,:,1) = [a22.*(eta - 1) - a21.*(xi - 1), a11.*(xi -1) - a12.*(eta - 1)]./[J, J];
Dphi(1:NT,:,2) = [a22.*(1 - eta) + a21.* xi,      -a11.*xi - a12.*(1 - eta)    ]./[J, J];
Dphi(1:NT,:,3) = [a22.*eta - a21.*xi,             a11.*xi - a12.*eta           ]./[J, J];
Dphi(1:NT,:,4) = [- a22.*eta - a21.*(1 - xi)    , a11.*(1 - xi) + a12.*eta     ]./[J, J];


