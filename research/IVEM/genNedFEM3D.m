function fem = genNedFEM3D(mesh,bcind,option)

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
    %bcind = [1,1,1,1,1,1];  
    ng = 4; nge = 2;
    option.ng = 4; option.nge = 2;
end
if nargin == 3
    ng = option.ng;
    nge = option.nge;
end
if ~isfield(option,'ng')
    option.ng = 4;
end

%% 1. Form fem p, t, and basis functions
fem = genNedbasFEM3D(mesh,1); option.ngauss = 4;

%% 2. Form Gaussian Quadrature
p = mesh.p; t = mesh.t;
X1 = p(t(:,1),:); X2 = p(t(:,2),:); X3 = p(t(:,3),:); X4 = p(t(:,4),:); 
gw = gaussWtetra(ng);
[gx,gy,gz] = gaussPtetra(X1,X2,X3,X4,ng);
A = tetraArea(X1,X2,X3,X4);
fem.gw = gw; fem.gx = gx; fem.gy = gy; fem.gz = gz; fem.area = A;
[gew,gex,gey,gez] = gaussPedge(p(mesh.e(:,1),:),p(mesh.e(:,2),:),nge);
fem.gew = gew; fem.gex = gex; fem.gey = gey; fem.gez = gez;

%% 3. Form bc mapper
[bc,mapper] = boundaryEdge3D(mesh.p,mesh.e,bcind);
fem.bc = bc; fem.mapper = mapper; fem.ng = option.ng;

%% 4. Form edge direction index
t = fem.t;
nt = size(t,1);
t_e_orit = zeros(nt,6);
e_ind = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
for i = 1:6
    d = mesh.p(t(:,e_ind(i,2)),:) - mesh.p(t(:,e_ind(i,1)),:);
    d_crect = mesh.p(mesh.e(mesh.t_e(:,i),2),:) - mesh.p(mesh.e(mesh.t_e(:,i),1),:);
    t_e_orit(:,i) = sign(sum(d.*d_crect,2));
end
fem.t_e_orit = t_e_orit;

%% 5. Assign the location of degrees of freedom
gdof = mesh.e;
fem.gdof = gdof;
g2ldof = mesh.t_e;
fem.g2ldof = g2ldof;
fem.bcind = bcind;




