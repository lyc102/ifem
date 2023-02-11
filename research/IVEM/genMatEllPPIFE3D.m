function matrix = genMatEllPPIFE3D(pde,mesh,fem,femI,femIF,PPtype,sig)
%% Generate global matrices and load vector of PPIFEM for 3D ellipic eq
%     -div(A grad u)  = f,    x\in \Omega
%      where A is a piecewise constant on Omega^+ and Omega^-.
% INPUTS:
% pde --- given data function from equation, e.g. 
%         pde.A --- diffusion coefficient
%         pde.f --- right hand side function
%         pde.gD --- Dirichlet boundary value function 
%         pde.one --- constant function 1.
% mesh --- mesh structure. 
% fem --- global degree of freedom of FEM 
% femI --- quadrature info on interface cells
% femIF --- quadrature info on interface faces, required in PPIFE.
% PPtype --- partial penalty type possible value 
%            'N' : Nonsymmetric PPIFE (e = 1)
%            'S' : Symmetric PPIFE (e = -1)
%            'I' : Incomplete PPIFE (e = 0)
% OUTPUTS:
% matrix.S --- stiffness matrix (w/o boundary condition)
% matrix.E --- consistence matrix (w/o boundary condition)
% matrix.P --- penalty matrix (w/o boundary condition)
% matrix.A --- final FEM matrix (after boundary condition)
% matrix.rhsF --- load vector (w/o boundary condition)
% matrix.rhsE --- consistence vector (=0 if no interface face on boundary)
% matrix.rhsP --- penalty vector (=0 if no interface face on boundary)
% matrix.f --- final RHS matrix (after boundary condition)

% Last Modified: 08/07/2020 by Xu Zhang

%% 1. Stiffness Matrix
S = globMatrixIFE3DStiff(pde.A,mesh,femI,fem,fem);
% This matrix generator is a little faster than doing it separately.
% Sx = globMatrixIFE3D(pde.A,mesh,femI,fem,[1,0,0],fem,[1,0,0]);
% Sy = globMatrixIFE3D(pde.A,mesh,femI,fem,[0,1,0],fem,[0,1,0]);
% Sz = globMatrixIFE3D(pde.A,mesh,femI,fem,[0,0,1],fem,[0,0,1]);
% S = Sx+Sy+Sz;

%% 2. Face Matrices
d1 = 0; j1 = 1; d2 = 1; j2 = 0;
E = globMatrixIFE3DFace(pde.one, pde.A, fem, fem, femIF, d1,j1, d2,j2);

d1 = 0; j1 = 1; d2 = 0; j2 = 1;
P = globMatrixIFE3DFace(pde.one, pde.one,fem,fem, femIF, d1,j1, d2,j2);

%% 3. Generate the Right Hand Side Vector
rhsF = globRHSIFE3D(pde.f, mesh, fem, femI, [0,0,0]);
rhsE = globRHSIFE3DFace(pde.A, pde.exactu, fem, femIF, 1);
rhsP = globRHSIFE3DFace(pde.one, pde.exactu, fem, femIF, 0);

%% 4. Dirichlet Boundary Conditions
if strcmp(PPtype,'N')
    e = 1; 
elseif strcmp(PPtype,'S')
    e = -1; 
elseif strcmp(PPtype,'I')
    e = 0;
end
h = mesh.p(2,1) - mesh.p(1,1);
Atotal = S - E + e*E' + (sig/h)*P;
rhs = rhsF + e*rhsE + (sig/h)*rhsP;
tu = feval(pde.gD,fem.p(:,1),fem.p(:,2),fem.p(:,3));
ub = tu;
ub(fem.mapper) = 0;
rhsB = Atotal*ub;
A = Atotal(fem.mapper,fem.mapper);
f = rhs(fem.mapper) - rhsB(fem.mapper);

%% 5. Outputs
matrix = struct('A',A,'f',f,'S',S,'E',E,'P',P,'rhsF',rhsF,'rhsE',rhsE,...
    'rhsP',rhsP);