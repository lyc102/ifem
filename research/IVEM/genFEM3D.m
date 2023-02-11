function fem = genFEM3D(mesh,femtype,bcind,option)

%% Usage: Form FEM global degrees of freedom on 3D tetrahedron mesh
%
% INPUTS:
% mesh --- a struct data contains rich mesh information.
% bcind --- [bc1,bc2,bc3,bc4,bc5,bc6] denotes the boundary condition of
%            bc1: x = xmin  bc2: x = xmax
%            bc3: y = ymin  bc4: y = ymax
%            bc5: z = zmin  bc6: z = zmax
%            if bc = 1: Dirichlet boundary condition
%               bc = 2: Neumann boundary condition
%               bc = 3: Robin boundary condition
% femtype --- the type of finite element methods
%             Conforming: P1, P2
%
% OUTPUTS:
% fem --- a struct data contains the following fields:
%         fem.p: (x,y,z) coordinate of each vertex w.r.t a global DoF
%         fem.t: indices of global DoF in each element
%         fem.bc: indices of global DoF on the boundary
%         fem.mapper: a vector to extract unknowns from global DoF.
%         fem.type: type of finite element methods
%         fem.ldof: number of local DoF on each element
%         fem.ng: number of Gaussian quadrature points on each element

% Last Modified: 08/07/2020 by Xu Zhang

%% 0. Inputs
if nargin == 2
    bcind = [1,1,1,1,1,1];  option.ng = 4;
end
if nargin == 3
    option.ng = 4;
end
if ~isfield(option,'ng')
    option.ng = 4;
end
ng = option.ng;

%% 1. Form fem p, t, and basis functions
if strcmp(femtype,'P1')
    fem = genP1FEM3D(mesh); option.ngauss = 4;    
elseif strcmp(femtype,'P2')
    fem = genP2FEM3D(mesh); option.ngauss = 11;
end

%% 2. Form Gaussian Quadrature
p = mesh.p; t = mesh.t;
X1 = p(t(:,1),:); X2 = p(t(:,2),:); X3 = p(t(:,3),:); X4 = p(t(:,4),:); 
gw = gaussWtetra(ng);
[gx,gy,gz] = gaussPtetra(X1,X2,X3,X4,ng);
A = tetraArea(X1,X2,X3,X4);
fem.gw = gw; fem.gx = gx; fem.gy = gy; fem.gz = gz; fem.area = A;

%% 3. Form bc mapper
[bc,mapper] = boundaryNode3D(mesh,fem.p,bcind);
fem.bc = bc; fem.mapper = mapper; fem.ng = option.ng; fem.bcind = bcind;
