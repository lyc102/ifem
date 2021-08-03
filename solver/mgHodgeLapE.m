function [x,info,Ai,Bi,BBi,Res,Pro] = mgHodgeLapE(A,b,node,elem,bdFlag,option,varargin)
%% MGHODGELAPE multigrid solver for the Schur complement of Hodge Laplacian
%
% [-Mv  G'] [p] = [f]
% [ G   C ] [u] = [g]
%
% The Schur complement A = C + G*DMvinv*G' is SPD. Multilevel coarsening
% is applied to the mesh and Galerkin projection used in the coarsen grids.
%
% Created by Jie Zhou on July 30, 2014.
%
% Reference: 
% L Chen, Y Wu, L Zhong, J Zhou. MultiGrid Preconditioners for Mixed Finite
% Element Methods of the Vector Laplacian. Journal of Scientific Computing
% 77 (1), 101-128.
%
% See also mgMaxwellsaddle, mgMaxwell
% 
% Revised by Long Chen on Feb 19, 2014, Aug 14, 2016.

t = cputime;
Ndof = size(b,1); % size of the system
N = max(elem(:)); % number of nodes
d = size(node,2); % dimension

%% Options
if ~exist('option','var') 
    option = []; 
end
if ~isfield(option,'smoothingstep')  % smoothing steps
    option.smoothingstep = 3;
end
if ~isfield(option,'smoothingratio')  % ratio of variable smoothing
    option.smoothingratio = 1.5;
end
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; 
N0 = option.N0; 
tol = option.tol;
maxIt = option.solvermaxit; 
mu = option.smoothingstep; 
smothingRatio = option.smoothingratio; % for variable smoothing
solver = option.solver; 
coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; 
setupflag = option.setupflag;

%% Set up multilevel structure
if setupflag == true
    %% Hierarchical meshes (uniform)
    HB = zeros(N,3);
    level = 20;
    NL(level+1) = N; % now NL(1:level) = 0;
    elemi = cell(level,1);
    nodei = cell(level,1);
    bdFlagi = cell(level,1);
    elemi{level} = elem;
    nodei{level} = node;
    bdFlagi{level} = bdFlag;
    for j = level: -1 : 2
        if d == 2
            [elemi{j-1},newHB,bdFlagi{j-1}] = uniformcoarsenred(elemi{j},bdFlagi{j});
        elseif d == 3
            [elemi{j-1},~,newHB,bdFlagi{j-1}] = uniformcoarsen3red(elemi{j},[],bdFlagi{j});
        end
        if ~isempty(newHB)            
            NL(j) = NL(j+1) - size(newHB,1); % update NL(k)
            HB(NL(j)+1:NL(j+1),1:3) = newHB(:,1:3);
            nodei{j-1} = nodei{j}(1:NL(j),:);
        else
        % no nodes are removed or it reaches the coarsest level
            NL = NL(j:end);       
            break; 
        end
        if size(elemi{j-1},1)< 2*N0
            NL = NL(j-1:end);
            break;
        end
    end
    level = length(NL)-1;    % actual level
    elemi = elemi(end-level+1:end);
    nodei = nodei(end-level+1:end);
    bdFlagi = bdFlagi(end-level+1:end);

    %% Matrices in each level
    % assemble Hodge Laplacian in each level
    Ai = cell(level,1);
    Ai{level} = A;    
    isFreeEdge = cell(level,1);
    isFreeEdge{level} = option.isFreeEdge;
    for j = level-1:-1:1
        if d == 2
            [Ai{j},isFreeEdge{j}] = HodgeLaplacianE_matrix(nodei{j},elemi{j},bdFlagi{j});
        elseif d == 3
            [Ai{j},isFreeEdge{j}] = HodgeLaplacian3E_matrix(nodei{j},elemi{j},bdFlagi{j});
        end
    end
    % transfer operator
    for j = level:-1:2
        if d == 2
            Pro{j-1} = transferedgered(elemi{j-1},elemi{j});
        elseif d == 3
            Pro{j-1} = transferedgered3(elemi{j-1},elemi{j});            
        end
        Pro{j-1}   = Pro{j-1}(isFreeEdge{j},isFreeEdge{j-1});
        Res{j} = Pro{j-1}';
        if isfield(option,'coarsematrix') && strcmp(option.coarsematrix,'G')
        % compute the coarse matrix by Galerkin projection
            Ai{j-1} = Res{j}*Ai{j}*Pro{j-1};
        end
    end    
    % smoother
    for j = level:-1:2
        switch option.smoother
            case 'GS'
                Bi{j} = tril(Ai{j});        % Forward Gauss-Seidel   B = D+L
                BBi{j} = triu(Ai{j});       % Backward Gauss-Seidel BB = D+U    
            case 'JAC'
                Bi{j} = spdiags(diag(Ai{j}),0,size(Ai{j},1),size(Ai{j},1));        %
                BBi{j} = Bi{j};             % Jacobi iteration    
        end
        if option.smoothingparameter~=1
            Bi{j} = Bi{j}/option.smoothingparameter;
            BBi{j} = BBi{j}/option.smoothingparameter;
        end
    end
end %% end of setup
if strcmp(solver,'NO')
    x = x0; flag = 0; itStep = 0; err = 0; time = cputime - t;
    info = struct('solverTime',time,'itStep',itStep,'error',err,...
                  'flag',flag,'stopErr',max(err(end,:)));        
    return
end

%% Multilevel structure is in the input
if setupflag == false
    Ai = varargin{1};
    Bi = varargin{2};
   BBi = varargin{3};
   Res = varargin{4};
   Pro = varargin{5};
   level = size(Ai,1);   
end

%% Multigrid cycles
if level == 1 
    % it is possible no multilevel structure is aviable since only edge
    % transfer for uniform coarsen for red refinement is implemented
    % then use AMG
    option.preconditioner = 'W';
    [x,info] = amg(A,b,option);  % then use amg
    return;
end
x = x0;
k = 1;
r = b-A*x;
nb = norm(b);
err = zeros(maxIt,2);
if nb > eps
    err(1,:) = norm(r)/nb;
else % b = 0
    err(1,:) = norm(r);
end
% solvers
switch solver
   case 'CG'
        if printlevel >= 2
          fprintf('Conjugate Gradient Method\n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)    
            % compute Br by MG
            Br = vcycle(r);
            % update tau, beta, and p
            rho = Br'*r;  % e'*ABA*e approximates e'*A*e
            if k == 1
                p = Br;
            else
                beta = rho/rho_old;
                p = Br + beta*p;
            end
            % update alpha, x, and r
            Ap = A*p;
            alpha = rho/(Ap'*p);
            r = r - alpha*Ap;
            x = x + alpha*p;
            rho_old = rho;
            k = k + 1;
            % compute err for the stopping criterion
        %     err(k,1) = alpha*sqrt(p'*Ap/(x'*A*x)); % increamental error in energy norm
            err(k,1) = sqrt(abs(rho/(x'*b))); % approximate relative error in energy norm
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u,  #nnz: %8.0u, MGCG iter: %2.0u, err = %8.4e\n',...
                         Ndof, nnz(A), k-1, max(err(k,:)));
            end
        end
    err = err(1:k,:);
    itStep = k-1;
    case 'VCYCLE'  
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
end

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
time = cputime - t;
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));    
if printlevel >= 1
    fprintf('Conjugate Multigrid for Hodge Laplacian \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(A), itStep, info.stopErr, time)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions: V-cycle and W-cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vcycle MG
    function Br = vcycle(r,J)      % solve equations Ae = r in each level  
        if nargin<=1
            J = level;
        end
        ri = cell(J,1);            % record residual in each level
        ei = cell(J,1);            % record err in each level
        ri{J} = r;
        for i = J:-1:2
            ei{i} = Bi{i}\ri{i};   % pre-smoothing
            for s = 1:ceil((smothingRatio^(level-i)))*mu-1   % variable smoothing
                ei{i} = ei{i} + Bi{i}\(ri{i}-Ai{i}*ei{i}); 
            end
            ri{i-1} = Res{i}*(ri{i} - Ai{i}*ei{i});

        end
        if strcmp(coarsegridsolver,'direct')
            ei{1} = Ai{1}\ri{1};    % direct solver in the coarest level
        else                        % iterative solver in the coarest level
            D = spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
            [ei{1}] = pcg(Ai{1},ri{1},1/size(Ai{1},1),1000,D);
        end
        for i = 2:J
            ei{i} = ei{i} + Pro{i-1}*ei{i-1};
            ei{i} = ei{i} + BBi{i}\(ri{i}-Ai{i}*ei{i});
            for s = 1:ceil((smothingRatio^(level-i)))*mu-1
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
        if strcmp(coarsegridsolver,'direct')
            e = Ai{1}\r;   % exact solver in the coaresest grid
        else                        % iterative solver in the coarest level
            D = spdiags(diag(Ai{1}),0,size(Ai{1},1),size(Ai{1},1));
            e = pcg(Ai{1},r,1/size(Ai{1},1),1000,D);
        end
        return
    end
    % fine grid pre-smoothing
    e = Bi{J}\r;        % pre-smoothing
    for s = 1:ceil((smothingRatio^(level-J)))*mu-1   % variable smoothing
        e = e + Bi{J}\(r-Ai{J}*e); 
    end
    % coarse grid correction twice
    rc = Res{J}*(r - Ai{J}*e);
    ec = wcycle(rc,J-1);
    ec = ec + wcycle(rc - Ai{J-1}*ec,J-1);
    e = e + Pro{J-1}*ec;
    % fine grid post-smoothing
    e = e + BBi{J}\(r-Ai{J}*e);
    for s = 1:ceil((smothingRatio^(level-J)))*mu-1
        e = e + BBi{J}\(r-Ai{J}*e); % post-smoothing
    end
    end
end