close all;

%%  sphere
surface = MolecularSurface('sphere.pqr');
cube = [-1, 1, -1, 1, -1, 1]*3;
h = 0.08;
[node,face,face2elem,interface] = interfacemesh3(cube, surface, h);
showmesh(node, interface.face);

%% two spheres
surface = MolecularSurface('twosphere.pqr');
cube = [-1, 1, -1, 1, -1, 1]*3;
h = 0.08;
[node,face,face2elem,interface] = interfacemesh3(cube, surface, h);
showmesh(node, interface.face);

%% ben molecular
surface = MolecularSurface('ben.pqr');
cube = [-1, 1, -1, 1, -1, 1]*3;
h = 0.08;
[node,face,face2elem,interface] = interfacemesh3(cube, surface, h);
showmesh(node, interface.face);