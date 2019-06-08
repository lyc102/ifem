function [x,info,Ai,Bi,BBi,Res,Pro,cl] = amg(A,b,option,varargin)
%% AMG algebraic multi-grid solvers
%
%   x = AMG(A,b) attempts to solve the system of linear equations A*x = b
%   for x using multigrid type solver. To acheive multigrid efficiency, a 
%   hierarchical 'grids' is generated from the graph of A. 
%
%   x = AMG(A,b,option) specifies options. In addition to options in <a href="matlab:help mg">mg</a>
%   the following two are added:
%   - option.theta: parameter used to define strong connectness
%   - option.coarsen: coarsening methods. Now the following is impelemented
%       * 'rs'  Ruge-Stuben coarsening
%       * 'c'   my modification of RS coarsening
%       * 'a'   aggregation based coarsening
%   - option.interpolation: interpolation methods used to construct prolongation. 
%       * 's' standard interpolation. Use the matrix A_fc as a weighted average
% of all connected coarse nodes.
%       * 't' two-points interpolation. Use at most two connected coarse
% nodes.
%       * 'a' aggegration interpolation. Use the strongest connected coarse
% node.
%       * 'sa' smoothing aggregation interpolation.
%
%   [x,info] = AMG(A,b) also returns information of the solver.
%   - info.flag:
%   	* 0: mg converged to the desired tolerance tol within maxIt iterations
%       * 1: mg iterated maxIt times but did not converge.
%   - info.itStep: the iteration number at which x was computed.
%   - info.time: the cpu time to get x
%   - info.err: the approximate relative error in the energy norm in
%   err(:,1) and the relative residual norm(b-A*x)/norm(b) in err(:,2). If
%   flag is 0, then max(err(end,:)) <= tol.
%   - info.stopErr: the error when iteration stops
%
%   Example:
%     load lakemesh
%     [node,elem] = uniformrefine(node,elem);
%     A = assemblematrix(node,elem);
%     [bdNode,bdEdge,isBdNode] = findboundary(elem);
%     A = A(~isBdNode,~isBdNode);
%     N = size(A,1);
%     b = ones(N,1)/N;
%     x = amg(A,b);
%
% See also mg
%
% Documentation in Help browser <a href="matlab:ifem amgdoc">ifem amgdoc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

t = cputime;
%% Initial setup of parameters
N = size(b,1);
Nb = size(b,2);               % number of bs. multiple rhs is allowed

%% Options
% Assign default values to unspecified parameters
if ~exist('option','var')
    option = []; 
end 
[theta,cmethod,interpolation,preconditioner,mu] = amgoption(option);
option = mgoptions(option,length(b));    % parameters
x0 = repmat(option.x0,1,Nb); 
N0 = option.N0; 
tol = option.tol; 
maxIt = option.solvermaxit; 
solver = option.solver; 
coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; 
setupflag = option.setupflag;
% additional parameter for amg
% Check for zero right hand side
if (norm(b) == 0)               % if rhs vector is all zeros
    x = zeros(N,1);             % then solution is all zeros
    flag = 0;  
    time = cputime - t;
    itStep = 0;                    
    err = 0;      
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));
    return
end

%% Minimal ordering for Gauss-Seidel smoothing
% So far no help on the performance
% ord = amd(A);  
% A = A(ord,ord);
% b = b(ord);

%% Set up multilevel structure
if setupflag == true
level = max(min(ceil(log2(N)/2-4),8),2);
Ai = cell(level,1);
Bi = cell(level,1);
BBi = cell(level,1);
Pro = cell(level,1);
Res = cell(level,1);
Ai{level} = A;
cl = 1; % coarsest level
for k = level:-1:2
    % coarsening
    switch upper(cmethod)
        case 'RS'
            [isC,As] = coarsenAMGrs(Ai{k},theta);
        case 'C'
            [isC,As] = coarsenAMGc(Ai{k},theta);
        case 'A' % aggregation
            [node2agg,As] = coarsenAMGa(Ai{k},theta);
%             N0 = 2.5*N0; % increase coarsen grid points for aggregation
    end
    % prolongation and restriction
    switch upper(interpolation)
        case 'S'  % standard
            [tempPro,tempRes] = interpolationAMGs(As,isC);
        case 'T'  % two points
            [tempPro,tempRes] = interpolationAMGt(As,isC);
        case 'A'  % one point
            [tempPro,tempRes] = interpolationAMGa(As,isC);
        case 'N'  % new by Xiaozhe
            [tempPro,tempRes] = interpolationAMGn(Ai{k},isC);
        case 'SA' % smooth aggregation
            [tempPro,tempRes] = interpolationAMGsa(As,node2agg);
    end
    % record operators
    Bi{k} = tril(Ai{k});   % smoother for pre-smoothing
    BBi{k} = triu(Ai{k});  % smoother for post-smoothing    
    Pro{k-1} = tempPro;
    Res{k} = tempRes;
    Ai{k-1} = tempRes*Ai{k}*tempPro;
    % check if reach the coarsest level
    if size(Ai{k-1},1)< N0 
        cl = k-1;
        break;
    end
end
end % end for setup

if strcmp(solver,'NO') % only set up the transfer matrices
    x = x0; flag = 0; itStep = 0; err = 0; time = cputime - t;
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));
    return
end

%% No need of set up
if setupflag == false
    Ai = varargin{1};
    Bi = varargin{2};
   BBi = varargin{3};
   Res = varargin{4};
   Pro = varargin{5};
   cl  = varargin{6};
   level = size(Ai,1);   
end

%% Krylov iterative methods use Multigrid-type Preconditioners
% set up preconditioner
switch preconditioner
    case 'V'
        prefunc = @vcycle;
        if printlevel >= 1
            fprintf('\n Algebraic Multigrid V-cycle Preconditioner with ')
        end
    case 'W'
        prefunc = @wcycle;
        if printlevel >= 1
            fprintf('\n Algebraic Multigrid W-cycle Preconditioner with ')
        end
    case 'F'
        prefunc = @fcycle;
        if printlevel >= 1
            fprintf('\n Algebraic Multigrid Full Cycle Preconditioner with ')
        end
    case 'BPX'
        % modify smoother
        for j = level:-1:2
            Di{j} = diag(Ai{j});
        end
        prefunc = @bpx;
        if printlevel >= 1
            fprintf('BPX Preconditioner with ')
        end
end
% initial set up
k = 1; 
x = x0;
r = b - A*x;
nb = max(sqrt(sum(b.^2,1)));
err = zeros(maxIt,2);
err(1,:) = max(sqrt(sum(r.^2,1)))/nb;
% solvers
switch solver
    case 'CG'
        if printlevel >= 1
            fprintf('Conjugate Gradient Method\n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)    
            % compute Br by MG
            Br = prefunc(r);
            % update tau, beta, and p
            rho = dot(Br,r);  % e'*ABA*e approximates e'*A*e
            if k == 1
                p = Br;
            else
                beta = rho./rho_old;
                p = Br + beta.*p;
            end
            % update alpha, x, and r
            Ap = A*p;
            alpha = rho./dot(Ap,p);
            r = r - alpha.*Ap;
            x = x + alpha.*p;
            rho_old = rho;
            k = k + 1;
            % compute err for the stopping criterion
        %     err(k,1) = alpha*sqrt(p'*Ap/(x'*A*x)); % increamental error in energy norm
            err(k,1) = max(sqrt(abs(rho./dot(x,b)))); % approximate relative error in energy norm
            err(k,2) = max(sqrt(sum(r.^2,1)))/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, MGCG iter: %2.0u, err = %8.4e\n',...
                         N, k-1, max(err(k,:)));
            end
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'VCYCLE'  
        if printlevel >= 1
            fprintf('\n Algebraic Multigrid Vcycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Vcylce MG
            Br = vcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = max(sqrt(abs(dot(Br,r)/dot(x,b)))); % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = max(sqrt(sum(r.^2,1)))/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, MG Vcycle iter: %2.0u, err = %8.4e\n',...
                         N, k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'WCYCLE'  
        if printlevel >= 1
            fprintf('\n Algebraic Multigrid Wcycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Vcylce MG
            Br = wcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = max(sqrt(abs(dot(Br,r)/dot(x,b)))); % approximate relative error in energy norm % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = max(sqrt(sum(r.^2,1)))/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, MG Wcycle iter: %2.0u, err = %8.4e\n',...
                         N, k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
    case 'FCYCLE'  
        if printlevel >= 1
            fprintf('\n Algebraic Multigrid Full Cycle Iteration \n')
        end
        while (max(err(k,:)) > tol) && (k <= maxIt)
            k = k + 1;
            % Step 2: Compute Br by one Full cycle MG
            Br = fcycle(r);
            % Step 3: Correct the solution
            x = x + Br;
            err(k,1) = max(sqrt(abs(dot(Br,r)/dot(x,b)))); % approximate relative error in energy norm % approximate relative error in energy norm
            % Step 1: Form residual r
            r = r - A*Br;
            err(k,2) = max(sqrt(sum(r.^2,1)))/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, MG Fcycle iter: %2.0u, err = %8.4e\n',...
                         N, k-1, max(err(k,:)));
            end            
        end
        err = err(1:k,:);
        itStep = k-1;
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
    case 'BICG'
        if printlevel >= 1
            fprintf('BiConjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = bicg(A,b,tol,maxIt,prefunc,[],x0);  
    case 'BICGSTAB'
        if printlevel >= 1
            fprintf('Stablilized BiConjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = bicgstab(A,b,tol,maxIt,prefunc,[],x0);  
    case 'BICGSTAB1'
        if printlevel >= 1
            fprintf('Stablilized BiConjugate Gradient Method\n')
        end
        [x,flag,err,itStep] = bicgstab1(A,b,tol,maxIt,prefunc,[],x0);  
end

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
time = cputime - t;
if printlevel >= 1
    fprintf('  nnz/N: %2.2f,   level: %2.0u,   coarse grid %2.0u,   nnz/Nc %2.2f\n',...
              nnz(A)/N, level-cl+1,size(Ai{cl},1),nnz(Ai{cl})/size(Ai{cl},1))
end
if printlevel >= 1
    fprintf('#dof: %7.0u,    iter: %2.0u,   err = %8.4e,   time = %4.3g s\n \n',...
                 N, itStep, max(err(end,:)), time)
end
if (flag == 1) && (printlevel > 0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));

%% Permute x back
% x(ord) = x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions vcycle, wcycle, fcycle, bpx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vcycle MG
    function Br = vcycle(r,J)  % solve equations Ae = r in each level  
    if nargin<=1
        J = level;
    end
    ri = cell(J,1);            % record residual in each level
    ei = cell(J,1);            % record err in each level
    ri{J} = r;
    for i = J:-1:cl+1
        ei{i} = Bi{i}\ri{i};   % pre-smoothing
        for s = 1:mu          % extra mu steps smoothing
            ei{i} = ei{i} + Bi{i}\(ri{i}-Ai{i}*ei{i}); 
        end
        ri{i-1} = Res{i}*(ri{i} - Ai{i}*ei{i});
    end
    if strcmp(coarsegridsolver,'direct')
        ei{cl} = Ai{cl}\ri{cl};    % direct solver in the coarest level
    else                        % iterative solver in the coarest level
        D = spdiags(diag(Ai{cl}),0,size(Ai{cl},1),size(Ai{cl},1));
        [ei{cl},flag] = pcg(Ai{cl},ri{cl},1/size(Ai{cl},1),1000,D);
    end
    for i = cl+1:J
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + BBi{i}\(ri{i}-Ai{i}*ei{i});
        for s = 1:mu
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
    if J == cl
        e = Ai{cl}\r;   % exact solver in the coaresest grid
        return
    end
    % fine grid pre-smoothing
    e = Bi{J}\r;   % pre-smoothing
    for s = 1:mu           % extra mu steps smoothing
        e = e + Bi{J}\(r-Ai{J}*e); 
    end
    rc = Res{J}*(r - Ai{J}*e);
    % coarse grid correction twice
    ec = wcycle(rc,J-1);
    ec = ec + wcycle(rc - Ai{J-1}*ec,J-1);
    % fine grid post-smoothing
    e = e + Pro{J-1}*ec;
    e = e + BBi{J}\(r-Ai{J}*e);
    for s = 1:mu
        e = e + BBi{J}\(r-Ai{J}*e); % post-smoothing
    end
    end

%% Fcycle MG
    function Br = fcycle(r)
    ri = cell(level,1);            % record residual in each level
    ei = cell(level,1);            % record err in each level
    ri{level} = r;
    for i = level:-1:cl+1
        ei{i} = vcycle(ri{i},i);   % pre-smoothing
        for s = 1:mu               % extra smoothing steps
            ei{i} = ei{i} + vcycle(ri{i}-Ai{i}*ei{i},i); % pre-smoothing
        end
        ri{i-1} = Res{i}*(ri{i} - Ai{i}*ei{i});
    end
    if strcmp(coarsegridsolver,'direct')
        ei{cl} = Ai{cl}\ri{cl};        % direct solver in the coarest level
    else                            % iterative solver in the coarest level
        D = spdiags(diag(Ai{cl}),0,size(Ai{cl},1),size(Ai{cl},1));
        [ei{cl},flag] = pcg(Ai{cl},ri{cl},1/size(Ai{cl},1),1000,D);
    end
    for i = cl+1:level
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + vcycle(ri{i}-Ai{i}*ei{i},i);
        for s = 1:mu   % extral smoothing
            ei{i} = ei{i} + vcycle(ri{i}-Ai{i}*ei{i},i); % post-smoothing
        end
    end
    Br = ei{level};    
    end

%% BPX preconditioner
    function Br = bpx(r)
    ri = cell(level,1);            % record residual in each level
    ei = cell(level,1);            % record err in each level
    % compute Br by BPX
    ri{level} = r;
    for i = level:-1:cl+1
        ei{i} = ri{i}./Di{i};       % Jacobi smoothing
        ri{i-1} = Res{i}*ri{i};     % restriction of the residual
    end
    if strcmp(coarsegridsolver,'direct')
        ei{1} = Ai{cl}\ri{1};         % direct solver in the coarest level
    else                           % iterative solver in the coarest level
        D = spdiags(diag(Ai{cl}),0,size(Ai{cl},1),size(Ai{cl},1));
        [ei{1},flag] = pcg(Ai{cl},ri{1},1/size(Ai{cl},1),1000,D);
    end
    for i = cl+1:level
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};  % prolongation of the correction
    end
    Br = ei{level};
    end

end