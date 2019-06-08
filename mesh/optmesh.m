function [node,elem] = optmesh(node,elem,rho)
%% OPTMESH optimize the input mesh
%
% [node,elem] = optmesh(node,elem) optimize the shape regularity of
% triangles in the input mesh (node,elem) and outputs a better mesh
% (node,elem). 
%
% [node,elem] = optmesh(node,elem,rho) accepts a non-uniform density given
% by an elementwise function rho. The default choice rho = 1/|t|. The
% quasi-uniform grids corresponds to rho=1. The density function rho can be
% given by a user specified function or a posteriori error estimator in the
% setting of adaptive finite element method.
%
% Example:
%  load airfoilperturbmesh
%  figure(1); subplot(1,2,1); 
%  showmesh(node,elem); title('original mesh');
%  figure(2); subplot(1,2,1); 
%  showmeshquality(node,elem); axis([0 1 0 2700]);
%  [node,elem] = optmesh(node,elem);
%  figure(1); subplot(1,2,2); 
%  showmesh(node,elem); title('smoothed mesh');
%  figure(2); subplot(1,2,2); 
%  showmeshquality(node,elem); axis([0 1 0 2700]);
%
% See also  bdsmoothing, meshsmoothing, edgeswap, rmisopoint, meshquality
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin<=2, rho = []; end
fprintf('Mesh quality before optimization \n')
meshquality(node,elem,rho);
node = meshsmoothing(node,elem,2,rho);
elem = edgeswap(node,elem);
[node,elem] = rmisopoint(node,elem);
node = meshsmoothing(node,elem,2,rho);
fprintf('Mesh quality after optimization \n')
meshquality(node,elem,rho);