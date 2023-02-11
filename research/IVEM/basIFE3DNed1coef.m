function [BAS1,BAS2] = basIFE3DNed1coef(vert,intpt,coef)

%% USAGE: Coefficients of P1 IFE Shape Function on a Tetrahedron
%
% INPUTS:
% vert --- 4-by-3 matrix stores the
%          vertices of the tetrahedron
%          [A1;A2;A3;A4] = [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4]
% intpt --- 3-by-3 or 4-by-3 matrix stores three intersection points
%          [D;E;F] = [xd yd zd; xe ye ze; xf yf zf];
%           the order is 1,2,3 or 2,3,4,5 of the index of edges
% coef --- [bt1,bt2] jump coefficient where bt1 is associated with A1.
%                                    bt2 is associated with A2 A3 A4
%
%     Type 1 Interface Element (3 intersection points), it means
%         A1 is on piece #1 whose coefficient is coef(1)
%         A2 A3 A4 are on piece #2 whose coefficient is coef(2)
%
%     Type 2 Interface Element (4 intersection points), it means
%         A1 A2 are on piece #1 whose coefficient is coef(1)
%         A3 A4 are on piece #2 whose coefficient is coef(2)
%
% OUTPUTS:
% BAS1 --- 6-by-6 matrix stores the coefficients of shape fun on piece #1
% BAS2 --- 6-by-6 matrix stores the coefficients of shape fun on piece #2
%          The ith basis: a * x + b on piece #j = 1 or 2
%              1st basis [a11,a12,a13,b11,b12,b13]
%              2nd basis [a21,a22,a23,b21,b22,b23]
%              3rd basis [a31,a32,a33,b31,b32,b33]
%              4th basis [a41,a42,a43,b41,b42,b43]
%              5th basis [a51,a52,a53,b51,b52,b53]
%              6th basis [a61,a62,a63,b61,b62,b63]
%
% Last Modified by Ruchi Guo on 19/11/20

%% 0. Initialization
bt1 = coef(1);  bt2 = coef(2); rb = bt2/bt1;
at1 = coef(3);  at2 = coef(4); ra = at2/at1;
 M = zeros(6,6); l = zeros(6,2); % the 2nd component is only associated with the shape fun on the piece2
if size(intpt,1) == 3
    D = intpt(1,:); E = intpt(2,:); F = intpt(3,:); C = 1/3*(D+E+F);
    M(1,1:3) = 1/2*(vert(1,:)+D); l(1,1) = norm(D-vert(1,:)); 
    M(1,4:6) = 1/2*(vert(2,:)+D); l(1,2) = norm(D-vert(2,:));
    M(2,1:3) = 1/2*(vert(1,:)+E); l(2,1) = norm(E-vert(1,:));
    M(2,4:6) = 1/2*(vert(3,:)+E); l(2,2) = norm(E-vert(3,:));
    M(3,1:3) = 1/2*(vert(1,:)+F); l(3,1) = norm(F-vert(1,:)); 
    M(3,4:6) = 1/2*(vert(4,:)+F); l(3,2) = norm(F-vert(4,:));
    M(4,4:6) = 1/2*(vert(2,:)+vert(3,:)); l(4,2) = norm(vert(2,:)-vert(3,:)); 
    M(5,4:6) = 1/2*(vert(2,:)+vert(4,:)); l(5,2) = norm(vert(2,:)-vert(4,:));
    M(6,4:6) = 1/2*(vert(3,:)+vert(4,:)); l(6,2) = norm(vert(3,:)-vert(4,:));
elseif size(intpt,1) == 4
    [intptU,~,~] = vert4to3(intpt);
    D = intptU(1,:); E = intptU(2,:); F = intptU(3,:); C = 1/3*(D+E+F);
    M(1,1:3) = 1/2*(vert(1,:)+vert(2,:)); l(1,1) = norm(vert(1,:)-vert(2,:));
    M(2,1:3) = 1/2*(vert(1,:)+intpt(1,:)); l(2,1) = norm(vert(1,:)-intpt(1,:));
    M(2,4:6) = 1/2*(vert(3,:)+intpt(1,:)); l(2,2) = norm(vert(3,:)-intpt(1,:));
    M(3,1:3) = 1/2*(vert(1,:)+intpt(2,:)); l(3,1) = norm(vert(1,:)-intpt(2,:));
    M(3,4:6) = 1/2*(vert(4,:)+intpt(2,:)); l(3,2) = norm(vert(4,:)-intpt(2,:));
    M(4,1:3) = 1/2*(vert(2,:)+intpt(3,:)); l(4,1) = norm(vert(2,:)-intpt(3,:));
    M(4,4:6) = 1/2*(vert(3,:)+intpt(3,:)); l(4,2) = norm(vert(3,:)-intpt(3,:));
    M(5,1:3) = 1/2*(vert(2,:)+intpt(4,:)); l(5,1) = norm(vert(2,:)-intpt(4,:));
    M(5,4:6) = 1/2*(vert(4,:)+intpt(4,:)); l(5,2) = norm(vert(4,:)-intpt(4,:));
    M(6,4:6) = 1/2*(vert(3,:)+vert(4,:)); l(6,2) = norm(vert(3,:)-vert(4,:));
end

t = zeros(6,6);
t(1,1:3) = vert(2,:)-vert(1,:); t(1,1:3) = t(1,1:3)/norm(t(1,1:3));
t(2,1:3) = vert(3,:)-vert(1,:); t(2,1:3) = t(2,1:3)/norm(t(2,1:3));
t(3,1:3) = vert(4,:)-vert(1,:); t(3,1:3) = t(3,1:3)/norm(t(3,1:3));
t(4,1:3) = vert(3,:)-vert(2,:); t(4,1:3) = t(4,1:3)/norm(t(4,1:3));
t(5,1:3) = vert(4,:)-vert(2,:); t(5,1:3) = t(5,1:3)/norm(t(5,1:3));
t(6,1:3) = vert(4,:)-vert(3,:); t(6,1:3) = t(6,1:3)/norm(t(6,1:3));

n = cross(E-D,F-D);
n = n/norm(n); % unit normal of plane DEF 
% L = @(x,y,z) dot([x-C(1),y-C(2),z-C(3)],n);

%% for the piece2
bas0 = zeros(6,6);
bas0(:,1) = (M(:,2).*t(:,3) - M(:,3).*t(:,2)).*l(:,1) + (M(:,5).*t(:,3) - M(:,6).*t(:,2)).*l(:,2);
bas0(:,2) = (M(:,3).*t(:,1) - M(:,1).*t(:,3)).*l(:,1) + (M(:,6).*t(:,1) - M(:,4).*t(:,3)).*l(:,2);
bas0(:,3) = (M(:,1).*t(:,2) - M(:,2).*t(:,1)).*l(:,1) + (M(:,4).*t(:,2) - M(:,5).*t(:,1)).*l(:,2);
bas0(:,4:6) = t(:,1:3).*((l(:,1)+l(:,2))*ones(1,3)); 

%% for the piece1
MM = M - [ones(6,1)*C, ones(6,1)*C];
bas1 = zeros(6,6);
bas1(:,1) = ((ra-1)*(MM(:,2).*t(:,3) - MM(:,3).*t(:,2)) +...
    (rb-1)*(t(:,1:3)*n')*(C(2)*n(3)-C(3)*n(2))).*l(:,1);
bas1(:,2) = ((ra-1)*(MM(:,3).*t(:,1) - MM(:,1).*t(:,3)) +...
    (rb-1)*(t(:,1:3)*n')*(C(3)*n(1)-C(1)*n(3))).*l(:,1);
bas1(:,3) = ((ra-1)*(MM(:,1).*t(:,2) - MM(:,2).*t(:,1)) +...
    (rb-1)*(t(:,1:3)*n')*(C(1)*n(2)-C(2)*n(1))).*l(:,1);
bas1(:,1:3) = bas1(:,1:3) -...
    (ra-1)*(sum(cross(ones(6,1)*n,MM(:,1:3)).*t(:,1:3),2).*l(:,1))*n;
bas1(:,4:6) = (rb-1)*(l(:,1).*(t(:,1:3)*n'))*n;


A = bas0 + bas1;

ss = A\eye(6);
ss = ss';
BAS2 = ss;
BAS1 = ss; 
BAS1(:,1:3) = ss(:,1:3) + (ra-1)*(ss(:,1:3)-(ss(:,1:3)*n')*n); 
BAS1(:,4:6) = ss(:,4:6) - (ra-1)*cross(ss(:,1:3)-(ss(:,1:3)*n')*n,ones(6,1)*C) +...
    (rb-1)*((cross(ss(:,1:3),ones(6,1)*C)+ss(:,4:6))*n')*n;
%BAS1 = BAS1'; BAS2 = BAS2';

return;
