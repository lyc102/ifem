close all;

%% circle
box = [ -1, 1, -1, 1];
h = 0.1;
phi = @(p) sum(p.^2, 2) - 0.5.^2;
[node,elem,interface] = interfacemesh(box,phi,h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');
pause

%% flower
[node,elem,interface] = interfacemesh(box,@phiflower,h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');
pause

%% heart
[node,elem,interface] = interfacemesh(box,@phiheart,0.5*h);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');