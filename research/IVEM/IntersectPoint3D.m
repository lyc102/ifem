function f = IntersectPoint3D(p1, p2, fun, option)

%% Usage: Use Bisection to find intersection points of a sequence of lines
%         and a surface
%
% INPUTS:
% p1 --- coordinates of first point(s).
% p2 --- coordinates of second point(s).
% fun --- function of the surface.
%
% OUTPUTS:
% f --- coordinates of the intersection point
%
% Last Modified: 06/29/2020 by Xu Zhang

%%
F1 = feval(fun, p1(:,1), p1(:,2), p1(:,3));
F2 = feval(fun, p2(:,1), p2(:,2), p2(:,3));
if isempty(sign(F1).*sign(F2) == 1)
    error('There is no intersection points!')
end

if size(nargin,1) == 3 && isfield(option,'surf') == 1
    
    for i = 1:60 %
        pM = (p1 + p2)/2; FM = feval(fun, pM(:,1), pM(:,2), pM(:,3));
        p1(FM.*F1 > 0,:) = pM(FM.*F1 > 0,:);
        p2(FM.*F2 > 0,:) = pM(FM.*F2 > 0,:);
        F1 = feval(fun, p1(:,1), p1(:,2), p1(:,3));
        F2 = feval(fun, p2(:,1), p2(:,2), p2(:,3));
    end
    f = (p1 + p2)/2;
    
else
    
    f = (-F2./(F1-F2)*ones(1,3)).*p1 + (F1./(F1-F2)*ones(1,3)).*p2;
    
end

% for i = 1:100 %
%     pM = (p1 + p2)/2; FM = feval(fun, pM(:,1), pM(:,2), pM(:,3));
%     p1(FM.*F1 > 0,:) = pM(FM.*F1 > 0,:);
%     p2(FM.*F2 > 0,:) = pM(FM.*F2 > 0,:);
%     F1 = feval(fun, p1(:,1), p1(:,2), p1(:,3));
%     F2 = feval(fun, p2(:,1), p2(:,2), p2(:,3));
% end
% f = (p1 + p2)/2;

return;
