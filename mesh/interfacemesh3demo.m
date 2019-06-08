close all;
% Get a molecular object 
% surface = MolecularSurface('adp.pqr')
% surface = MolecularSurface('twosphere.pqr');
surface = MolecularSurface('ben.pqr');
cube = [-1, 1, -1, 1, -1, 1]*3;
h = 0.08;
[node, elem, interfaceData] = interfacemesh3(cube, surface, h);
showmesh(node, interfaceData.interface);
disp([cube, h]);
disp(surface)
