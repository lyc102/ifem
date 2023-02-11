function [p,t] = genMesh3DRectPT(domain, nx, ny, nz)

%% Usage: nodes and elements info of a uniform cubic mesh
%
% INPUTS:
% domain --- cubic domain = [xmin, xmax, ymin, ymax, zmin, zmax].
% nx --- the number of uniform partition in x direction.
% ny --- the number of uniform partition in y direction.
% nz --- the number of uniform partition in z direction.

% OUTPUTS:
% p --- np-by-3 vector: (x,y,z) coordinates of each node.
% t --- nt-by-8 vector: eight-node indices for each cube, ordered as follows
%                    8----7
%                   /|   /|
%                  5-|--6 |
%                  | |  | |
%                  | 4--|-3 
%                  |/   |/
%                  1----2
% Last Modified: 08/07/2020 by Xu Zhang

%% 0. Initial Setting
xmin = domain(1); xmax = domain(2); hx = (xmax-xmin)/nx;
ymin = domain(3); ymax = domain(4); hy = (ymax-ymin)/ny;
zmin = domain(5); zmax = domain(6); hz = (zmax-zmin)/nz;
x = xmin:hx:xmax; nnx = nx+1;
y = ymin:hy:ymax; nny = ny+1;
z = zmin:hz:zmax; nnz = nz+1;

%% 1. Form p
X = repmat(x',nny*nnz,1);
Y = repmat(reshape(repmat(y,nnx,1),nnx*nny,1),nnz,1);
Z = reshape(repmat(z,nnx*nny,1),nnx*nny*nnz,1);
p = [X,Y,Z];

%% 2. Form t
% 2.1 form a 2D rectangular mesh first
tmp = reshape((1:nnx*ny)',nnx,ny);
tmp(end,:) = [];
c = reshape(tmp,nx*ny,1); 
t2D = [c,c+1,c+nx+2,c+nx+1];

% 2.2 Form the 3D partition
t = zeros(nx*ny*nz,8);
for k=1:nz
    t(1+(k-1)*(nx*ny):k*nx*ny,:) = [t2D+(k-1)*nnx*nny, t2D+k*nnx*nny];
end
