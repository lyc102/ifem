
function [c] = tricircumcenter3d(node,elem)

% http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html

v1 = node(elem(:,1),:);
v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:);

ba = v2 - v1;
ca = v3 - v1;


baLength = sum(ba.^2,2);
caLength = sum(ca.^2,2);


crossbc = cross(ba,ca,2);

d = 0.5 ./sum(crossbc.^2,2);


xc = ((baLength.*ca(:,2) - caLength.*ba(:,2)).*crossbc(:,3) - (baLength.*ca(:,3) - caLength.*ba(:,3)).*crossbc(:,2)).*d;
yc = ((baLength.*ca(:,3) - caLength.*ba(:,3)).*crossbc(:,1) - (baLength.*ca(:,1) - caLength.*ba(:,1)).*crossbc(:,3)).*d;
zc = ((baLength.*ca(:,1) - caLength.*ba(:,1)).*crossbc(:,2) - (baLength.*ca(:,2) - caLength.*ba(:,2)).*crossbc(:,1)).*d;

c = [xc, yc, zc];

c = v1 + c;


