function [BAS1,BAS2] = basIFE3DP1coef(vert,intpt,coef)

%% USAGE: Coefficients of P1 IFE Shape Function on a Tetrahedron
%
% INPUTS:
% vert --- 4-by-3 matrix stores the
%          vertices of the tetrahedron
%          [A1;A2;A3;A4] = [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4]
% intpt --- 3-by-3 matrix stores three intersection points
%          [D;E;F] = [xd yd zd; xe ye ze; xf yf zf];
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
% BAS1 --- 4-by-4 matrix stores the coefficients of shape fun on piece #1
% BAS2 --- 4-by-4 matrix stores the coefficients of shape fun on piece #2
%          The ith basis: ci1 + ci2 x + ci3 y + ci4 z on piece #j = 1 or 2
%              1st basis [c11,c12,c13,c14]
%              2nd basis [c21,c22,c23,c24]
%              3rd basis [c31,c32,c33,c34]
%              4th basis [c41,c42,c43,c44]
%
% Last Modified by Xu Zhang on 08/07/20

%% 0. Initialization
bt1 = coef(1);  bt2 = coef(2); r = bt1/bt2;
if size(intpt,1) == 3
    D = intpt(1,:); E = intpt(2,:); F = intpt(3,:); C = 1/3*(D+E+F);
elseif size(intpt,1) == 4
    [intptU,~,~] = vert4to3(intpt);
    D = intptU(1,:); E = intptU(2,:); F = intptU(3,:); C = 1/3*(D+E+F);
else
    stp = 1;
end

n = cross(E-D,F-D);
n = n/norm(n); % unit normal of plane DEF 
% L = @(x,y,z) dot([x-C(1),y-C(2),z-C(3)],n);

%% 1. Impose the jump conditions:
bas1 = eye(4);
bas2 = bas1 + (r-1)*[0,n]'.*repmat([dot(-C,n),n],4,1);

%% 2. Impose the nodal conditions:
P = [ones(4,1),vert];
if size(intpt,1) == 3
    A = [P(1,:)*bas1'; P(2,:)*bas2'; P(3,:)*bas2'; P(4,:)*bas2'];
elseif size(intpt,1) == 4
    A = [P(1,:)*bas1'; P(2,:)*bas1'; P(3,:)*bas2'; P(4,:)*bas2'];
end
ss = A\eye(4);
BAS1 = ss'*bas1;
BAS2 = ss'*bas2;
