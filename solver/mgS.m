function [x,info] = mgS(A,b,elem,option,varargin)
%% MGS multigrid-type solvers
%
%   x = MGS(A,b,elem) attempts to solve the system of linear equations A*x =
%   b for x using multigrid type solver. To acheive multigrid efficiency,
%   the mesh must be of bisection type. Inside mg, an coarsening algorithm
%   on bisection grids is applied. See <a href="matlab:ifem coarsendoc">doc coarsen</a> for the coarsening algorithm. 
% 
%   It is a simplified version of |mg| such that it can serve as a template
%   for coding multigrid methods using iFEM.
%
% Example:
%
%   [node,elem] = squaremesh([0 1 0 1],0.5);
%   for k = 1:8
%     [node,elem] = uniformrefine(node,elem);
%   end
%   pde.f = inline('p(:,1).*p(:,2)','p');
%   pde.g_D = inline('zeros(size(p,1),1)','p');
%   option.solver = 'notsolve';
%   [u,Du,eqn] = Poisson(node,elem,pde,[],option);
%   fprintf('\n Number of unknowns: %8.0u\n',length(eqn.b))
%   tic; display('Direct solver'); u = eqn.A\eqn.b; toc;
%   tic; x = mgS(eqn.A,eqn.b,elem); toc;
%   format shorte
%   fprintf('Difference between direct and mg solvers %0.2g \n',norm(u-x));
%
% See also mg
%
% Documentation in Help browser <a href="matlab:ifem mgdoc">ifem mgdoc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tic;

%% Size of systems
Ndof = length(b);                  % number of dof
N = max(elem(:));                  % number of nodes
% NT = size(elem,1);                 % number of elements
dim = size(elem,2)-1;

%% Options
% Assign default values to unspecified parameters
if ~exist('option','var'), 
    option = []; 
end
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; 
N0 = option.N0; 
tol = option.tol;
maxIt = option.solvermaxit; 
mu = option.smoothingstep; 
solver = option.solver; 
preconditioner = option.preconditioner;
coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; 
% setupflag = option.setupflag;

%% Fixed dof
isFixDof = [];
if isfield(option,'freeDof') % freeDof is given
   isFreeDof = false(N,1);   
   isFreeDof(option.freeDof) = true;
   NA = size(A,1);
   if NA > length(option.freeDof)
       isFixDof = true(NA,1);
       isFixDof(isFreeDof) = false;
   end
else % Find free dof and eliminate isolated dof
    deg = sum(spones(A));  % degree of the matrix graph
    isFreeDof = false(Ndof,1);
    isFreeDof(deg>1) = true;
    isFixDof = ~isFreeDof; % degree = 1 is isolated vertex
    isFixDof(deg == 0) = false; % degree = 0 the matrix is singular
end
if any(isFixDof) % a bigger matrix is given
    xD = zeros(Ndof,1);
    xD(isFixDof) = b(isFixDof)./diag(A(isFixDof,isFixDof));
    A = A(isFreeDof,isFreeDof);% shrink to freedof only
end
if length(b) > sum(isFreeDof)  % a bigger b is given
    b = b(isFreeDof);          % shrink to freedof only
end
if length(x0) > sum(isFreeDof) % a bigger x0 is given
    x0 = x0(isFreeDof);        % shrink to freedof only
    option.x0 = x0;
end
isFreeNode = isFreeDof(1:N); % free nodes (for P1)

%% Hierarchical Structure of Mesh
%  hierarchical structure for linear element.
if dim == 2  % 2-D
    [HB, NL, level] = HBstructure(elem,N0); 
end
if dim == 3  % 3-D
    if nargin > 4
        HBmesh = varargin{end}; % need HB for bisection refinement
        [HB, NL, level] = HBstructure3(elem,N0,HBmesh);
    else % no HBmesh is given. only work for the red uniform refinement
        [HB, NL, level] = HBstructure3(elem,N0);
    end
end

%% Transfer operators between multilevel meshes for P1 element
% standard prolongation and restriction operator for P1 element
[Pro,Res] = transferoperator(HB,NL,isFreeNode); 
if Ndof > N
    if ~exist('auxPro','var')
        disp('The current element is not supported by mg');
    else
        Pro{level} = auxPro;
        Res{level+1} = Pro{level}'; 
        level = level + 1; % add one more level from P1 to current element
    end
end
clear HB auxPro

%% Matrices in each level
Ai = cell(level,1);
Ai{level} = A;    
for j = level:-1:2
    Ai{j-1} = Res{j}*Ai{j}*Pro{j-1};           % Ac = Res*Af*Pro
    switch option.smoother
        case 'GS'
            Bi{j} = tril(Ai{j});        % Forward Gauss-Seidel   B = D+L
            BBi{j} = triu(Ai{j});       % Backward Gauss-Seidel BB = D+U    
        case 'JAC'
            Bi{j} = spdiags(diag(Ai{j}),0,size(Ai{j},1),size(Ai{j},1));        %
            BBi{j} = Bi{j};       % Jacobi iteration    
    end
    if option.smoothingparameter~=1
        Bi{j} = Bi{j}/option.smoothingparameter;
        BBi{j} = BBi{j}/option.smoothingparameter;
    end
end

%% MG cycles
% initial set up
k = 1; 
x = x0;
r = b - A*x;
nb = norm(b);
err = zeros(maxIt,2);
if nb > eps  % nb is non-zero
    err(1,:) = norm(r)/nb; 
else
    err(1,:) = norm(r);
end
switch solver
    case 'VCYCLE'  
        if printlevel >= 1
            fprintf('Multigrid Vcycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Vcylce MG
            Br = vcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,1) = sqrt(abs(Br'*r/(x'*b))); % approximate relative error in energy norm
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u,  #nnz: %8.0u, MG Vcycle iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'WCYCLE'  
        if printlevel >= 1
            fprintf('Multigrid Wcycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Vcylce MG
            Br = wcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = sqrt(abs(Br'*r/(x'*b))); % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u,  #nnz: %8.0u, MG Wcycle iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
end

%% Krylov space method
% set up preconditioner
switch preconditioner
    case 'V'
        prefunc = @vcycle;
        if printlevel >= 1
            fprintf('Multigrid V-cycle Preconditioner with ')
        end
    case 'W'
        prefunc = @wcycle;
        if printlevel >= 1
            fprintf('Multigrid W-cycle Preconditioner with ')
        end
end
switch solver
    case 'CG'
        if printlevel >= 1
            fprintf('Conjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = pcg(A,b,tol,maxIt,prefunc,[],x0);          
    case 'MINRES'
        if printlevel >= 1
            fprintf('Minimum Residual Method \n')
        end
        [x,flag,err,itStep] = minres(A,b,tol,maxIt,prefunc,[],x0);  
    case 'GMRES'
        if printlevel >= 1
            fprintf('General Minimum Residual Method\n')
        end
        if isfield(option,'restart')
            restart = option.restart;
        else
            restart = min(N,10);
        end
        [x,flag,err,itStep] = gmres(A,b,restart,tol,maxIt,prefunc,[],x0);
        itStep = (itStep(1)-1)*restart + itStep(2);
end

%% Modify x to include fix dof
if any(isFixDof)
    xD(isFreeDof) = x;
    x = xD;
end

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
time = toc;
if printlevel >= 2
    fprintf('#dof: %8.0u, level: %2.0u,   coarse grid %2.0u, #nnz: %8.0u\n',...
              Ndof, level, size(Ai{1},1), nnz(Ai{1}))
end
if printlevel >= 1
    fprintf('#dof: %8.0u,  #nnz: %8.0u, smoothing: %2.0u, iter: %2.0u,   err = %8.4e,   time = %4.2g s\n',...
                 Ndof, nnz(A), mu, itStep, max(err(end,:)), time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions vcycle, wcycle, fcycle, bpx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vcycle MG
    function Br = vcycle(r,J)        % solve equations Ae = r in each level  
    if nargin<=1
        J = level;
    end
    ri = cell(J,1);            % record residual in each level
    ei = cell(J,1);            % record err in each level
    ri{J} = r;
    for i = J:-1:2
        ei{i} = Bi{i}\ri{i};   % pre-smoothing
        for s = 1:mu-1           % extra mu-1 steps smoothing
            ei{i} = ei{i} + Bi{i}\(ri{i}-Ai{i}*ei{i}); 
        end
        ri{i-1} = Res{i}*(ri{i} - Ai{i}*ei{i});
    end
    if strcmp(coarsegridsolver,'direct')
        ei{1} = Ai{1}\ri{1}; % direct solver in the coarest level
    else                         % iterative solver in the coarest level
        D = spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
        [ei{1},flag] = pcg(Ai{1},ri{1},1/size(Ai{1},1),1000,D);
    end
    for i = 2:J
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + BBi{i}\(ri{i}-Ai{i}*ei{i});
        for s = 1:mu-1
            ei{i} = ei{i} + BBi{i}\(ri{i}-Ai{i}*ei{i}); % post-smoothing
        end
    end
    Br = ei{J};
    end

%% Wcycle MG
    function e = wcycle(r,J)        % solve equations Ae = r in each level  
    if nargin<=1
        J = level;
    end
    if J == 1
        e = Ai{J}\r;   % exact solver in the coaresest grid
        return
    end
    % fine grid pre-smoothing
    e = Bi{J}\r;        % pre-smoothing
    for s = 1:mu-1        % extra mu-1 steps smoothing
        e = e + Bi{J}\(r-Ai{J}*e); 
    end
    % coarse grid correction twice
    rc = Res{J}*(r - Ai{J}*e);
    ec = wcycle(rc,J-1);
    ec = ec + wcycle(rc - Ai{J-1}*ec,J-1);
    e = e + Pro{J-1}*ec;
    % fine grid post-smoothing
    e = e + BBi{J}\(r-Ai{J}*e);
    for s = 1:mu-1
        e = e + BBi{J}\(r-Ai{J}*e); % post-smoothing
    end
    end

end