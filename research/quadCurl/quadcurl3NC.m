function [soln,eqn] = quadcurl3NC(node,elem,bdFlag,pde,option,HB)
%% QUADCURL3NC: nonconforming approximation to the quad curl problem in 3D
%
% [u,eqn,info] = quadcurl3NC(node,elem,bdFlag,pde,option)
% computes the ND0-(CR-P0)-ND0 nonconforming approximations.
%
%  curl curl w = f
%  div w = 0
%
%  -\Delta phi + \nabla p = curl w
%  div phi = 0
%
%  curl curl u = curl phi
%  div u = g
%
% with Dirichlet boundary condition on \partial \Omega
% w x n = 0
% phi = 0
% u x n = 0
% 
% phi = curl u 
% (curl u) x n = 0 
% (curl u ) \dot n  = rot_{\Gamma} u = 0 by a density argument
% 
% see also Maxwell, Stokes3CRP0
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Set up optional input arguments
if ~exist('option','var'), option = []; end

%% Sort elem to ascend ordering
% [elem,bdFlag] = sortelem3(elem,bdFlag);
% [elem2face,face] = dof3face(elem);
% [elem2edge,edge] = dof3edge(elem);
% NF = size(face,1); NT = size(elem,1); NE = size(edge,1);

%% first curl problem
pdeCurl1.J = pde.quadcurlu;
pdeCurl1.g_D = @(p) zeros(size(p,1),3);
% option.solver = 'mg';
wh = Maxwellsaddle(node,elem,bdFlag,pdeCurl1,option,HB);

%% RHS of Stokes
[curlwh, volume, curlPhi] = curlu3(node,elem,wh);

%% Stokes part
pdeStokes.f = curlwh;
% pdeStokes.g_D = @(p) zeros(size(p,1),3);
pdeStokes.g_D = pde.curlu;
% option.solver = 'diag';
soln = Stokes3CRP0(node,elem,bdFlag,pdeStokes,option,HB);

%% last Maxwell equation's rhs
[elem2face,face] = dof3face(elem);
NF = size(face,1); 
NT = size(elem,1);
phih = reshape(soln.u,NF,3);
phihCenter = zeros(NT,3); % at each elem center
for j = 1:3 % each component
    phihj = phih(:,j);
    phihj2elem = phihj(elem2face);
    phihCenter(:,j) = sum(phihj2elem,2)/4;
    % Crouzeix-Raviart shape function 1-3\lambda = 1/4 at center
end

[elem2edge,edge] = dof3edge(sort(elem,2)); % sort for edge element    
NE = size(edge,1);
bt = zeros(NT,6);
for j = 1:6
    bt(:,j) = dot(curlPhi(:,:,j),phihCenter,2).*volume;
end
fCurl = accumarray(elem2edge(:),bt(:),[NE 1]);

%%
pdeCurl2.J = fCurl;
% pdeCurl2.g_D = @(p) zeros(size(p,1),3);
pdeCurl2.g_D = pde.exactu;
pdeCurl2.g = pde.g;
% option.solver = 'mg';
[maxsoln, eqnCurl] = Maxwellsaddle(node,elem,bdFlag,pdeCurl2,option,HB);

%% outputs
soln = struct('u',maxsoln.u,'w',wh, 'phi', phih, 'curlw', curlwh);
eqn = struct('A',eqnCurl.A,'G',eqnCurl.G);