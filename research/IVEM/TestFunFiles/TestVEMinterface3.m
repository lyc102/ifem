surface = MolecularSurface('twosphere.pqr');
cube = [-1, 1, -1, 1, -1, 1];
h = 2/10;
tic
[node,elem,face,face2elem,interface] = interfacemesh3(cube, surface, h);
toc
%showmesh(node, interface.face);
pde = twosphereinterfacedata(pi/3, 1, 10);
interface.face = face;
interface.face2elem = face2elem;
tic
[u, w, AE, AI, info] = interfacePoisson3VEM(node, elem, pde, interface, []);
toc