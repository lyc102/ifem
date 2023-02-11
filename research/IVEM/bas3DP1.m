function bas = bas3DP1(p,t)

%% USAGE: coefficient of Lagrange P1 basis for tetrahedral mesh
%
% INPUTS:
% p --- np-by-3 vector
% t --- nt-by-4 vector
%
% OUTPUTS:
% bas --- nt-4-4 matrix stores the coefficient of basis value(s)
%         The k-i-j th entry denotes the j-th basis function on k-th cell
%         with the i-th coefficient where c_i (i=1,2,3,4) stands for
%                 c_1 + c_2*x + c_3*y + c_4*z

% Last Modified: 08/07/2020 by Xu Zhang

%%
x1 = p(t(:,1),1); y1 = p(t(:,1),2); z1 = p(t(:,1),3); 
x2 = p(t(:,2),1); y2 = p(t(:,2),2); z2 = p(t(:,2),3); 
x3 = p(t(:,3),1); y3 = p(t(:,3),2); z3 = p(t(:,3),3); 
x4 = p(t(:,4),1); y4 = p(t(:,4),2); z4 = p(t(:,4),3); 

A = -x2.*y3.*z1 + x2.*y4.*z1 + x1.*y3.*z2 - x1.*y4.*z2 + x2.*y1.*z3 ...
    -x1.*y2.*z3 + x1.*y4.*z3 - x2.*y4.*z3 + x4.*(-y2.*z1 + y3.*z1 + ...
    y1.*z2 - y3.*z2 - y1.*z3 + y2.*z3) + (-x2.*y1 + x1.*y2 - x1.*y3 + ...
    x2.*y3).*z4 + x3.*(-y4.*z1 - y1.*z2 + y4.*z2 + y2.*(z1 - z4) + y1.*z4);

nt = size(t,1);  bas = zeros(nt,4,4);
bas(:,1,1) = (-x4.*y3.*z2 + x3.*y4.*z2 + x4.*y2.*z3 - x2.*y4.*z3 - x3.*y2.*z4 + x2.*y3.*z4);
bas(:,2,1) =-(-y3.*z2 + y4.*z2 + y2.*z3 - y4.*z3 - y2.*z4 + y3.*z4);
bas(:,3,1) = (-x3.*z2 + x4.*z2 + x2.*z3 - x4.*z3 - x2.*z4 + x3.*z4);
bas(:,4,1) =-(-x3.*y2 + x4.*y2 + x2.*y3 - x4.*y3 - x2.*y4 + x3.*y4);
bas(:,1,2) =-(-x4.*y3.*z1 + x3.*y4.*z1 + x4.*y1.*z3 - x1.*y4.*z3 - x3.*y1.*z4 + x1.*y3.*z4);
bas(:,2,2) = (-y3.*z1 + y4.*z1 + y1.*z3 - y4.*z3 - y1.*z4 + y3.*z4);
bas(:,3,2) =-(-x3.*z1 + x4.*z1 + x1.*z3 - x4.*z3 - x1.*z4 + x3.*z4);
bas(:,4,2) = (-x3.*y1 + x4.*y1 + x1.*y3 - x4.*y3 - x1.*y4 + x3.*y4);
bas(:,1,3) = (-x4.*y2.*z1 + x2.*y4.*z1 + x4.*y1.*z2 - x1.*y4.*z2 - x2.*y1.*z4 + x1.*y2.*z4);
bas(:,2,3) =-(-y2.*z1 + y4.*z1 + y1.*z2 - y4.*z2 - y1.*z4 + y2.*z4);
bas(:,3,3) = (-x2.*z1 + x4.*z1 + x1.*z2 - x4.*z2 - x1.*z4 + x2.*z4);
bas(:,4,3) =-(-x2.*y1 + x4.*y1 + x1.*y2 - x4.*y2 - x1.*y4 + x2.*y4);
bas(:,1,4) =-(-x3.*y2.*z1 + x2.*y3.*z1 + x3.*y1.*z2 - x1.*y3.*z2 - x2.*y1.*z3 + x1.*y2.*z3);
bas(:,2,4) = (-y2.*z1 + y3.*z1 + y1.*z2 - y3.*z2 - y1.*z3 + y2.*z3);
bas(:,3,4) =-(-x2.*z1 + x3.*z1 + x1.*z2 - x3.*z2 - x1.*z3 + x2.*z3);
bas(:,4,4) = (-x2.*y1 + x3.*y1 + x1.*y2 - x3.*y2 - x1.*y3 + x2.*y3);
bas = bas./A;