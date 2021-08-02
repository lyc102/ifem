% To Do: add more examples

close all;

%% circle
box = [ -1, 1, -1, 1];
h = 0.1;
phi = @(p) sum(p.^2, 2) - 0.5.^2;
[node,elem,interface] = interfacemesh(box,phi,h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');

%% flower
[node,elem,interface] = interfacemesh(box,@phiflower,h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');

%%
[node,elem,interface] = interfacemeshdoc(box,@phi4,h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');

%%
function [z,n] = phiflower(p) % level set function
petal = 5;
xc = 0;
yc = 0;
r0 = 0.5;
r1 = 0.2;
theta = atan2((p(:,2) - yc),(p(:,1) - xc));
r = r0 + r1*sin(petal*theta);

z = (p(:,2)-yc).^2 + (p(:,1)-xc).^2 - r.^2;

if nargout == 2
    rp = petal*r1*cos(petal*theta);
    n1 = rp.*cos(theta) - r.*sin(theta);
    n2 = rp.*sin(theta) + r.*cos(theta);
    l = sqrt(n1.^2+n2.^2);
    n = [n2./l, -n1./l];
end
end

function z = phi4(p)
x = p(:,1) - 0.02; y = p(:,2) -0.02;
z = (x.^2 + y.^2 - 0.38).^3 - x.^2.*y.^3;
end

function z = phi5(p)
x = p(:,1); y = p(:,2);
r = x.^2 + y.^2;
theta = atan2(y,x);
isNeg = theta < 0;
theta(isNeg) = theta(isNeg) + 2*pi;

x1 = 16*sin(theta).^3;
y1 = 13*cos(theta) - 5*cos(2*theta) - 2*cos(3*theta) - cos(4*theta);
z = r - (x1.^2 + y1.^2);
end