
domain = [-1,1,-1,1,-1,1];
% bm = 1; bp = 100;
% rx = 1; ry = 0.075; rz = 3;
% pde = elli3DorthocircIntf(bm,bp,rx,ry,rz);
%intf = @(x,y,z) x.^2+y.^2+z.^2-1;epsm = 8.85*10^(-3)*2;
epsp = 8.85*10^(-3);
sigm = 1;
sigp = 0.01;
mum = (4*pi)*3;
mup = (4*pi);
x0=0; y0=0; z0=-0.3; r1=0.2; r2=pi/5; omega = 1; a = omega*sqrt(epsp*mup); b= 150; intPt = -1;
pde = TorusTimeInitial3(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt);

fimplicit3(pde.intf,domain)
material metal

colormap(jet)

%%%%%%%%%%%%%%%%%%%%%%%
nx = 10;  h=(domain(2) - domain(1))/nx;
ny = nx;
nz = nx;

mesh = genMesh3D(domain, nx, ny, nz);
mesh = enrichMesh3D(mesh,2); % Mesh detail level = 1 (for IFE).
mesh = genIntfMesh3D(mesh,pde.intf);
meshI = genIVmesh(mesh);

mdpt = (mesh.p(mesh.t(:,1),:) + mesh.p(mesh.t(:,2),:) +...
    mesh.p(mesh.t(:,3),:) + mesh.p(mesh.t(:,4),:))/4;
tidshow = find(mdpt(:,2)>0);
tetramesh(mesh.t(tidshow,:),mesh.p,'FaceAlpha',0)
hold on
trisurf(meshI.iface,meshI.node(:,1),meshI.node(:,2),meshI.node(:,3)) %'edgecolor','none'
box on
%axis off
axis equal
set(gca,'XTick',[],'YTick',[],'ZTick',[])

view(50,30)