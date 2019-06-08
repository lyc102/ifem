% SOLVER
%
% Files
%   adgs                    - augmented Distributive Gauss-Seidel method for solving Stokes equations
%   amg                     - algebraic multi-grid solvers
%   amgoption               - AMGOPTIONS default options for algebraic multigrid solver
%   coarsenAMGc             - coarsen the graph of A.
%   coarsenAMGrs            - load lakemesh % For debug
%   gradmatrix              - matrix for the gradient of a nodal linear element
%   HBstructure             - reconstruct a hierachical structure of the mesh.
%   HBstructure3            - reconstructs a hierachical structure of a 3-D mesh.
%   interpolationAMGn       - INTERPOLATIONAMGS construct prolongation and restriction matrices
%   interpolationAMGs       - construct prolongation and restriction matrices
%   mg                      - multigrid-type solvers
%   mgMaxwell               - multigrid-type solver for Maxwell equations.
%   mgnew                   - MG multigrid-type solvers
%   mgoptions               - default options for multigrid solver
%   mgStokesP1P0Twolevel    - MGSTOKES multigrid-type solvers
%   mgStokesP1P0TwolevelNew - MGSTOKES multigrid-type solvers
%   node2edgematrix         - matrix to transfer nodal element to the lowest order edge element
%   node2edgematrix1        - matrix to transfer nodal element to the linear edge element
%   PFGMRES                 - Preconditioned Flexible General Minimal Residual Method (Right Preconditioner)
%   StokesDGSSmooth         - STOKESDGSSmooth Matrix-vecotr form DGS relaxation for Stokes eqns
%   transfermg              - solves saddle point problem discretized from mixed FEM.
%   transferoperator        - transfer operators between multilevel meshes
%   uniformcoarsen          - uniform coarsening
%   uniformcoarsen3         - uniform coarsening in 3-D
%   uniformcoarsen3red      - UNIFORMCOARSENRED uniform coarsening of red refinement
%   uniformcoarsenred       - uniform coarsening of red refinement
%   uzawa                   - augmented Uzawa method for solving Stokes equations
%   uzawapcg                - solves saddle point problem discretized from mixed FEM.
