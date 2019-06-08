function [u,p,edge,A,eqn,info] = Stokes(node,elem,pde,bdFlag,option)
%% STOKES Stokes equation
%
%   [u,p] = STOKES(node,elem,pde,bdFlag) use continuous and piecewise
%   quadratic element for velocity u = [u1,u2] and continuous and piceswise
%   linear element for pressure p to approximate the Stokes equations
% 
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%   with 
%       Dirichlet boundary condition        u = g_D  on \Gamma_D, 
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.
%
%   The mesh is given by node and elem and the boundary edge is given by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, and mu.
%
%  [u,p,edge,A,eqn,info] = STOKES(node,elem,pde,bdFlag,option) also
%  returns the edge structure of the mesh, the matrix A for the Laplace
%  operator, the eqn strcuture, and the info structure for solver information. 
%     - eqn.AD: modified Laplace matrix
%     - eqn.BD: modified negative divergence matrix
%     - eqn.f:  vector f
%     - eqn.g:  vector g
%
%  The solution [u,p] satisfy the saddle point problem:
%     [AD BD'][u] = [f]
%     [BD  0 ][p] = [g]
%
%  In the modified matrices AD, BD, the Dirichlet boundary condition of u
%  is built-into the equation as I_{\Gamma_D} u_D = g_D. The Neumann
%  boundary condition is built-into f by boundary integrals.
%
%  Users can chose other elements and specifies the solver options.
%     - option.fem: various elements for Stokes equations
%       * 'P2P1'    P2-P1 Taylor-Hood elements
%       * 'P2P0'    P2-P0 elements
%       * 'isoP2P0' P1(refined mesh)-P0 elements
%       * 'isoP2P1' P1(refined mesh)-P1 elements
%       * 'CRP0'    nonconforming P1 (CR)-P0 elements
%       * 'Mini'    P1+B-P1 elements
%
%     - option.solver: various solvers for Stokes system
%       * 'direct'  the built in direct solver \ (mldivide)
%       * 'uzawa'   inexact (augmented) Uzawa's method  
%       * 'mg'      multigrid based on DGS smoother
%       * 'blkdiag' block diagonal preconditioned MINRES
%       * 'blktri'  block triangular preconditioned GMRES
%   The default setting is to use the direct solver for small size problems
%   and 'blkdiag' mehthod for large size problems.
%
% Example
%   squareStokes;
%
% See also Poisson, StokesP2P0, StokesisoP2P0
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option.fem = 'P2P1'; end
if ~exist('bdFlag','var'), bdFlag = []; end
if ~isfield(option,'solver'), option.solver = 'blkdiag'; end

switch upper(option.elemType)
    case 'P2P1'
        [u,p,edge,A,eqn,info] = StokesP2P1(node,elem,pde,bdFlag,option);
    case 'P2P0'
        [u,p,edge,A,eqn,info] = StokesP2P0(node,elem,pde,bdFlag,option);
    case 'ISOP2P1'
        [u,p,edge,A,eqn,info] = StokesisoP2P1(node,elem,pde,bdFlag,option);
    case 'ISOP2P0'
        [u,p,edge,A,eqn,info] = StokesisoP2P0(node,elem,pde,bdFlag,option);
    case 'CRP0'
        [u,p,edge,A,eqn,info] = StokesCRP0(node,elem,pde,bdFlag,option);
    case 'CRP1'
        [u,p,edge,A,eqn,info] = StokesCRP1(node,elem,pde,bdFlag,option);
    case 'MINI'
        [u,p,edge,A,eqn,info] = StokesMini(node,elem,pde,bdFlag,option);
    case 'P1BP1'
        [u,p,edge,A,eqn,info] = StokesP1bP1(node,elem,pde,bdFlag,option);
    case 'RTP0'
        [u,p,w,edge,A,eqn,info] = StokesRT0(node,elem,pde,bdFlag,option);        
    case 'BDMP0'
        [u,p,w,edge,A,eqn,info] = StokesBDM1B(node,elem,pde,bdFlag,option);        
end