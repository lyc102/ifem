function matrix = genMatCurlIFE3D(pde,mesh,fem,femI)
%% Generate global matrices and load vector of FEM for 3D ellipic eq
%     -div(A grad u)  = f,    x\in \Omega
% INPUTS:
% pde --- given data function from equation, e.g. 
%         pde.A --- diffusion coefficient
%         pde.f --- right hand side function
%         pde.gD --- Dirichlet boundary value function 
%         pde.one --- constant function 1.
% mesh --- mesh structure. 
% fem --- global degree of freedom of FEM 
%
% OUTPUTS:
% matrix.S --- stiffness matrix (w/o boundary condition)
% matrix.A --- final FEM matrix (after boundary condition)
% matrix.rhsF --- load vector (w/o boundary condition)
% matrix.f --- final RHS matrix (after boundary condition)

% Last Modified: 08/07/2020 by Xu Zhang

%% 1. Stiffness Matrix
S = globMatrixNedIFE3D(pde.A,1,mesh,femI,fem,fem);
M = globMatrixNedIFE3D(pde.B,0,mesh,femI,fem,fem);

%% 2. Generate the Right Hand Side Vector
rhsF1 = globNedRHSIFE3D(pde.f1, mesh, fem, femI, 0, 1);
rhsF2 = globNedRHSIFE3D(pde.f2, mesh, fem, femI, 0, 2);
rhsF3 = globNedRHSIFE3D(pde.f3, mesh, fem, femI, 0, 3);
rhsF = rhsF1 + rhsF2 + rhsF3;

%% 3. Dirichlet Boundary Conditions
Atotal = S + M;
% non-interface edges
tu1 = sum(feval(pde.exactu1,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tu2 = sum(feval(pde.exactu2,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tu3 = sum(feval(pde.exactu3,fem.gex,fem.gey,fem.gez).*fem.gew,2);
% interface edges
eID = find(mesh.eLoc<0);
tu1I1 = sum(feval(pde.um1,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
tu2I1 = sum(feval(pde.um2,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
tu3I1 = sum(feval(pde.um3,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
tu1I2 = sum(feval(pde.up1,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
tu2I2 = sum(feval(pde.up2,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
tu3I2 = sum(feval(pde.up3,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
tu1(eID) = tu1I1 + tu1I2; tu2(eID) = tu2I1 + tu2I2; tu3(eID) = tu3I1 + tu3I2; 

tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
tgt = tgt./sum(tgt.^2,2).^(1/2);
tu = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
ub = tu;
ub(fem.mapper) = 0;
rhsB = Atotal*ub;
A = Atotal(fem.mapper,fem.mapper);
f = rhsF(fem.mapper) - rhsB(fem.mapper);

%% Outputs
matrix = struct('A', A, 'f', f, 'S', S, 'rhsF',rhsF,'tu',tu);