function [phi, Dphi, J] = hexbasis(node,elem,GaussPt)
%% HEXBASIS the basis and their derivatives and the det of Jacobi.
% 
%   [phi, Dphi, J] = hexbasis(node,elem,GaussPts) compute the values of 
%   bilinear basis and their derivatives and Jacobi det at a gauss point GaussPt in a
%   hex.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

NT = size(elem,1);

xi = GaussPt(1);
eta = GaussPt(2);
theta = GaussPt(3);

%% basis function and their derivatives
% for a point (x, y) in real element,  there exist a point in reference
% element (xi, eta), satisfy
%    x = x_1*phi_1 + x_2*phi_2 + x_3*phi_3 + x_4 *phi_4 +
%        x_5*phi_5 +x_6*phi_6 + x_7*phi_7 + x_8 *phi_8
%    y = y_1*phi_1 + y_2*phi_2 + y_3*phi_3 + y_4 *phi_4+
%        y_5*phi_5 +y_6*phi_6 + y_7*phi_7 + y_8 *phi_8
%    z = z_1*phi_1 + z_2*phi_2 + z_3*phi_3 + z_4 *phi_4+
%        z_5*phi_5 +z_6*phi_6 + z_7*phi_7 + z_8 *phi_8


g1 = [xi, 1 - xi];
g2 = [eta, 1 - eta];
g3 = [theta, 1 - theta];

phi = zeros(8,1);
phi(1) = g1(2)*g2(2)*g3(2);
phi(2) = g1(1)*g2(2)*g3(2);
phi(3) = g1(1)*g2(1)*g3(2);
phi(4) = g1(2)*g2(1)*g3(2);
phi(5) = g1(2)*g2(2)*g3(1);
phi(6) = g1(1)*g2(2)*g3(1);
phi(7) = g1(1)*g2(1)*g3(1);
phi(8) = g1(2)*g2(1)*g3(1);

rDphi = zeros(8,3);
rDphi(1,:) = [ -g2(2)*g3(2), -g1(2)*g3(2), -g1(2)*g2(2)];
rDphi(2,:) = [  g2(2)*g3(2), -g1(1)*g3(2), -g1(1)*g2(2)];
rDphi(3,:) = [  g2(1)*g3(2),  g1(1)*g3(2), -g1(1)*g2(1)];
rDphi(4,:) = [ -g2(1)*g3(2),  g1(2)*g3(2), -g1(2)*g2(1)];
rDphi(5,:) = [ -g2(2)*g3(1), -g1(2)*g3(1),  g1(2)*g2(2)];
rDphi(6,:) = [  g2(2)*g3(1), -g1(1)*g3(1),  g1(1)*g2(2)];
rDphi(7,:) = [  g2(1)*g3(1),  g1(1)*g3(1),  g1(1)*g2(1)];
rDphi(8,:) = [ -g2(1)*g3(1),  g1(2)*g3(1),  g1(2)*g2(1)];

X = node(:,1);
Y = node(:,2);
Z = node(:,3);

JM = zeros(NT,3,3);
JM(:,:,1) = X(elem)*rDphi;
JM(:,:,2) = Y(elem)*rDphi;
JM(:,:,3) = Z(elem)*rDphi;

J = JM(:,1,1).*JM(:,2,2).*JM(:,3,3) + JM(:,2,1).*JM(:,3,2).*JM(:,1,3)+ ...
    JM(:,3,1).*JM(:,1,2).*JM(:,2,3) - JM(:,3,1).*JM(:,2,2).*JM(:,1,3)- ...
    JM(:,1,1).*JM(:,3,2).*JM(:,2,3) - JM(:,2,1).*JM(:,1,2).*JM(:,3,3);

InvJM = zeros(NT,3,3);
InvJM(:,:,1) = [JM(:,2,2).*JM(:,3,3) - JM(:,3,2).*JM(:,2,3), JM(:,3,1).*JM(:,2,3)-JM(:,2,1).*JM(:,3,3), JM(:,2,1).*JM(:,3,2) - JM(:,3,1).*JM(:,2,2)];
InvJM(:,:,2) = [JM(:,3,2).*JM(:,1,3) - JM(:,1,2).*JM(:,3,3), JM(:,1,1).*JM(:,3,3)-JM(:,3,1).*JM(:,1,3), JM(:,3,1).*JM(:,1,2) - JM(:,1,1).*JM(:,3,2)];
InvJM(:,:,3) = [JM(:,1,2).*JM(:,2,3) - JM(:,2,2).*JM(:,1,3), JM(:,2,1).*JM(:,1,3)-JM(:,1,1).*JM(:,2,3), JM(:,1,1).*JM(:,2,2) - JM(:,2,1).*JM(:,1,2)];

InvJM(:,:,1) = InvJM(:,:,1)./[J,J,J];
InvJM(:,:,2) = InvJM(:,:,2)./[J,J,J];
InvJM(:,:,3) = InvJM(:,:,3)./[J,J,J];

Dphi = zeros(NT, 3, 8);

for i = 1:8
  Dphi(1:NT,:,i) = [InvJM(:,:,1)*rDphi(i,:)', InvJM(:,:,2)*rDphi(i,:)', InvJM(:,:,3)*rDphi(i,:)'];
end

