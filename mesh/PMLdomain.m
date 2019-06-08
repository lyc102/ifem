function [node,elem] = PMLdomain(node,h)
%Give a domain including node and  elem. We can construct a suitable
%pmldomain.
x0 = min(node(:,1));
x1 = max(node(:,1));
y0 = min(node(:,2));
y1 = max(node(:,2));

Lx = (x1 - x0)/10;% The thickness of PML along X;
Ly = (y1 - y0)/10;      % The thickness of PML along Y;

[node,elem]  = squaremesh([x0-Lx,x1+Lx,y0-Ly,y1+Ly],h);






















