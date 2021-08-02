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

%% heart
[node,elem,interface] = interfacemeshdoc(box,@phiheart,0.2*h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');