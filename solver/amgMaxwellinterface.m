function [x,info] = amgMaxwellinterface(A,b,node,edge,option)
%% AMGMAXWELL algebraic multigrid solver for Maxwell equations.
% 
% x = amgMaxwell(A,b,node,edge) attempts to solve the system of
% linear equations Ax = b using multigrid type solver. The linear system is
% obtained by either the first or second family of linear edge element
% discretization of the Maxwell equation; See <a href="matlab:ifem
% coarsendoc">doc Maxwell</a>.
%
% amgMaxwell is more algebraic than mgMaxwell but still requires geometric
% information node and edge. Grapha Laplacian of vertices are used as
% auxiliary Poisson operator and amg is used as Poisson solver.
%
% Input 
%   -  A: curl(alpha curl) + beta I
%   -  b: right hand side
%   -  node,edge: mesh information
%   -  options: extra structures
%
% By default, the HX preconditioned PCG is used which works well for
% symmetric positive definite matrices (e.g. arising from eddy current
% simulation). For symmetric indefinite matrices (e.g. arising from time
% harmonic Maxwell equation), set option.solver = 'minres' or 'bicg' or
% 'gmres' to try other Krylov method with HX preconditioner.
%
% See also mg, mgMaxwell, mgMaxwellsaddle
%
% Reference: 
% R. Hiptmair and J. Xu, Nodal Auxiliary Space Preconditioning in
% H(curl) and H(div) Spaces. SIAM J. Numer. Anal., 45(6):2483-2509, 2007.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

t = cputime;
%% Initial check
Ndof = length(b);                  % number of dof
N = size(node,1);                  % number of nodes
NE = size(edge,1);                 % number of edge;
dim = size(node,2);      
% Assign default values to unspecified parameters
if ~exist('option','var')
    option = []; 
end
option = mgoptions(option,length(b));    % parameters
x0 = option.x0; 
% N0 = option.N0; 
tol = option.tol; 
% mu = option.smoothingstep; preconditioner = option.preconditioner; coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; %setupflag = option.setupflag;
maxIt = 400;   % increase the default step (200) for Maxwell equations
% tol = 1e-7;    % reset the tol
% Check for zero right hand side
if (norm(b) == 0)                  % if rhs vector is all zeros
    x = zeros(Ndof,1);             % then solution is all zeros
    flag = 0;                      
    itStep = 0;                    
    err = 0;  
    time = cputime - t;
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));        
    return
end

%% Transfer operators from nodal element to edge element
if isfield(option,'isBdEdge')
    isBdEdge = option.isBdEdge;
else
    deg = sum(spones(A(1:NE,1:NE)),2);
    isBdEdge = (deg == 1);
end
if Ndof == NE        % lowest order edge element
    II = node2edgematrix(node,edge,isBdEdge);
elseif Ndof >= 2*NE  % first or second order edge element
    II = node2edgematrix1(node,edge,isBdEdge);
end
IIt = II';
[grad,isBdNode] = gradmatrix(edge,isBdEdge);
gradt = grad';

%% Free node 
% isFreeNode = false(N,1);
% isFreeNode(edge(isBdEdge,:)) = true;
% isFreeNode = true(N,1);

%% Auxiliary Poisson matrix
%   -  A: curl(alpha curl) + beta I
%   - AP: - div(alpha grad) + beta I
%   - BP: - div(beta grad)

if isfield(option,'AP')
    AP = option.AP;
else
    % build graph Laplacian to approximate AP
    edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
    edgeLength = sqrt(sum(edgeVec.^2,2));
    if isfield(option,'alpha') % resacle by the magnetic coefficients
        if isreal(option.alpha) && (length(option.alpha) == NE)
           alpha = option.alpha;  
        else % option.alpha is a function 
           edgeMiddle = (node(edge(:,2),:) + node(edge(:,1),:))/2; 
           alpha = option.alpha(edgeMiddle);         
        end
    else
        alpha = 1;
    end
    % edge weight: h*alpha
    % AP = gradt*spdiags(edgeLength.*alpha,0,NE,NE)*grad;
    i1 = (1:NE)'; j1 = double(edge(:,1)); s1 = sqrt(edgeLength.*alpha); 
    i2 = (1:NE)'; j2 = double(edge(:,2)); s2 = -s1;
    isFreeEdge = ~isBdEdge;
    G = sparse([i1(isFreeEdge);i2(isFreeEdge)],...
               [j1(isFreeEdge);j2(isFreeEdge)],...
               [s1(isFreeEdge);s2(isFreeEdge)],NE,N);
    AP = G'*G;
    % lumped mass matrix: h^3*beta
    if isfield(option,'beta') % resacle by the dielectric coefficients
        if isreal(option.beta) && (length(option.beta) == NE)
           beta = option.beta;  
        else % option.beta is a function 
           edgeMiddle = (node(edge(:,2),:) + node(edge(:,1),:))/2; 
           beta = option.beta(edgeMiddle);         
        end
    else
        beta = 1;
    end
    M = accumarray(edge(:),repmat((edgeLength.^3).*beta,2,1),[N 1]);
    AP = AP + spdiags(M,0,N,N);
end
% BP is Galerkin projection to the free node space
% boundary nodes
bdidx = zeros(N,1); 
if isempty(isBdNode)  % Neumann boundary condition
   bdidx = 1e-6;
else
   bdidx(isBdNode) = 1;
end
Tbd = spdiags(bdidx,0,N,N);
BP = gradt*A(1:NE,1:NE)*grad + Tbd;
% BP = gradt*spdiags(edgeLength.*beta,0,NE,NE)*grad + Tbd;

%% Transfer operators between multilevel meshes
setupOption.solver = 'NO';
[x,info,APi,Ri,RRi,ResAP,ProAP] = amg(AP,ones(N,1),setupOption); %#ok<ASGLU>
[x,info,BPi,Si,SSi,ResBP,ProBP] = amg(BP,ones(N,1),setupOption); %#ok<ASGLU>
D = diag(A);

%% Krylov iterative methods with HX preconditioner
k = 1;
err = 1;
switch upper(option.outsolver)
    case 'CG'
        if printlevel>=1
            fprintf('Conjugate Gradient Method using HX preconditioner \n');
        end
        x = x0;
        r = b - A*x;
        nb = norm(b);
        err = zeros(maxIt,2);
        err(1,:) = norm(r)/nb; 
        while (max(err(k,:)) > tol) && (k <= maxIt)
            % compute Br by HX preconditioner
            Br = HXpreconditioner(r); 
            % update tau, beta, and p
            rho = r'*Br;  % r'*Br = e'*ABA*e approximates e'*A*e
            if k==1
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
            % compute err for the stopping criterion
            k = k + 1;
            err(k,1) = sqrt(abs(rho/(x'*b))); % approximate relative error in energy norm
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, HXCG iter: %2.0u, err = %12.8g\n',...
                         Ndof, k, max(err(k,:)));
            end
        end
        err = err(1:k,:);
        itStep = k-1;
        if k > maxIt || (max(err(end,:))>tol)
            flag = 1;
        else
            flag = 0;
        end
    case 'MINRES'
        fprintf('Minimum Residual Method with HX preconditioner \n')
        [x,flag,err,itStep] = minres(A,b,tol,maxIt,@HXpreconditioner,[],x0);         
%         x = minres(A,b,tol,maxIt,@HXpreconditioner,[],x0);         
    case 'GMRES'
        fprintf('General Minimum Residual Method with HX preconditioner \n')
        [x,flag,err,itStep] = gmres(A,b,10,tol,maxIt,@HXpreconditioner,[],x0);
        itStep = (itStep(1) -1)*10 + itStep(2);
end

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('#dof: %8.0u,   #nnz: %8.0u,   iter: %2.0u,   err = %8.4e,   time = %4.2g s\n',...
                 Ndof, nnz(A), itStep, max(err(end,:)), time)
end
if  (flag == 1) && (printlevel>0)
    fprintf('NOTE: the iterative method does not converge');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions HXpreconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Br = HXpreconditioner(r)
    %% 1. Smoothing in the finest grid of the original system
%     eh = triu(A)\(D.*(tril(A)\r));  % Gauss-Seidal. less iteration steps
    eh = 0.75*r./D;  % Jacobi method. less computational time

    %% 2. Correction in the auxiliary spaces
    % Part1: II*(AP)^{-1}*II^t
    rc = reshape(IIt*r(1:size(IIt,2)),N,dim);   % transfer to the nodal linear element space
    level = size(APi,1);
    ri = cell(level,1);        
    ei = cell(level,1);        
    ri{level} = rc;            
    for i = level:-1:2
        ei{i} = Ri{i}\ri{i};   
        ri{i-1} = ResAP{i}*(ri{i}-APi{i}*ei{i});
    end
    ei{1} = APi{1}\ri{1};      
    for i = 2:level
        ei{i} = ei{i} + ProAP{i-1}*ei{i-1};
        ei{i} = ei{i} + RRi{i}\(ri{i} - APi{i}*ei{i});
    end
    eaux = [II*reshape(ei{level},dim*N,1); zeros(Ndof-size(IIt,2),1)]; % transfer back to the edge element space
    % Part2:  grad*(BP)^{-1}*grad^t
    level = size(BPi,1);
    ri{level} = gradt*r(1:NE);          % transfer the residual to the null space    
    for i = level:-1:2
        ei{i} = Si{i}\ri{i};   
        ri{i-1} = ResBP{i}*(ri{i}-BPi{i}*ei{i});
    end
    ei{1} = BPi{1}\ri{1};      
    for i=2:level
        ei{i} = ei{i} + ProBP{i-1}*ei{i-1};
        ei{i} = ei{i} + SSi{i}\(ri{i}-BPi{i}*ei{i});
    end
    eaux = eaux + [grad*ei{level}; zeros(Ndof-NE,1)]; % transfer back to the edge element space
    % combined correction in the finest grid
    Br = eh + eaux;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
