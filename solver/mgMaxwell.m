function [x,flag,itStep,err,time] = mgMaxwell(A,b,AP,BP,node,elem,edge,HBmesh,isBdEdge,option)
%% MGMAXWELL multigrid-type solver for Maxwell equations.
% 
% x = mgMaxwell(A,b,AP,BP,node,elem,edge,HBmesh,isBdEdge) attempts to solve
% the system of linear equations Ax = b using multigrid type solver. The
% linear system is obtained by either the first or second family of linear
% edge element discretization of the Maxwell equation; See <a href="matlab:ifem coarsendoc">doc Maxwell</a>.
%
% Input 
%   -  A: curl(mu^(-1)curl) + epsilon I
%   -  b: right hand side
%   - AP: - div(mu^{-1}grad) + |epsilon| I
%   - BP: - div(|epsilon|grad)
%   - node,elem,edge,HBmesh,isBdEdge: mesh information
%   - options: type help mg
%
% By default, the HX preconditioned PCG is used which works well for
% symmetric positive definite matrices (e.g. arising from eddy current
% simulation). For symmetric indefinite matrices (e.g. arising from time
% harmonic Maxwell equation), set option.solver = 'minres' or 'bicg' or
% 'gmres' to try other Krylov method with HX preconditioner.
%
% See also mg, mgMaxwellsaddle
%
% Reference: 
% R. Hiptmair and J. Xu, Nodal Auxiliary Space Preconditioning in
% H(curl) and H(div) Spaces. SIAM J. Numer. Anal., 45(6):2483-2509, 2007.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tic;
%% Initial check
Ndof = length(b);                  % number of dof
N = max(elem(:));                  % number of nodes
NE = size(edge,1);                 % number of edge;
dim = size(elem,2)-1;
% Assign default values to unspecified parameters
if nargin < 10 
    option = []; 
end 
option = mgoptions(option,length(b));    % parameters
x0 = option.x0; N0 = option.N0; tol = option.tol; 
solver = option.solver; 
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
    time = toc;
    return
end

%% Hierarchical Structure of Mesh
%  hierarchical structure for linear element.
if dim == 2  % 2-D
   [HB, NL, level] = HBstructure(elem,N0);        
end
if dim == 3  % 3-D
   [HB, NL, level] = HBstructure3(elem,N0,HBmesh);
end

%% Transfer operators from nodal element to edge element
if Ndof == NE        % lowest order edge element
    II = node2edgematrix(node,edge,isBdEdge);
elseif Ndof >= 2*NE  % first or second order edge element
    II = node2edgematrix1(node,edge,isBdEdge);
end
IIt = II';
grad = gradmatrix(edge,isBdEdge);
gradt = grad';

%% Free node 
% isFreeNode = false(N,1);
% isFreeNode(edge(isBdEdge,:)) = true;
isFreeNode = true(N,1);

%% Transfer operators between multilevel meshes
[Pro,Res] = transferoperator(HB,NL,isFreeNode);
APi = cell(level,1);
BPi = cell(level,1);
APi{level} = AP;    
BPi{level} = BP;                
for j = level:-1:2
    APi{j-1} = Res{j}*APi{j}*Pro{j-1};           % Ac = Res*Af*Pro    
    BPi{j-1} = Res{j}*BPi{j}*Pro{j-1};  
    Ri{j} = tril(APi{j});     % smoother for AP
    RRi{j} = triu(APi{j});   
    Si{j} = tril(BPi{j});     % smoother for BP   
    SSi{j} = triu(BPi{j});       
end
D = diag(A);                   
clear HB

%% Krylov iterative methods with HX preconditioner
k = 1;
err = 1;
switch upper(solver)
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
if k > maxIt || (max(err(end,:))>tol)
    flag = 1;
else
    flag = 0;
end
time = toc;
if printlevel >= 1
    fprintf('#dof: %8.0u,   #nnz: %8.0u,   iter: %2.0u,   err = %8.4e,   time = %4.2g s\n',...
                 Ndof, nnz(A), itStep, max(err(end,:)), time)
end
if  (flag == 1) && (printlevel>0)
    fprintf('NOTE: the iterative method does not converge');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions HXpreconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Br = HXpreconditioner(r)
    %% 1. Smoothing in the finest grid of the original system
    eh = triu(A)\(D.*(tril(A)\r));  % Gauss-Seidal. less iteration steps
%     eh = r./D;  % Jacobi method. less computational time

    %% 2. Correction in the auxiliary spaces
    % Part1: II*(AP)^{-1}*II^t
    rc = reshape(IIt*r(1:size(IIt,2)),N,dim);   % transfer to the nodal linear element space
    ri = cell(level,1);        
    ei = cell(level,1);        
    ri{level} = rc;            
    for i = level:-1:2
        ei{i} = Ri{i}\ri{i};   
        ri{i-1} = Res{i}*(ri{i}-APi{i}*ei{i});
    end
    ei{1} = APi{1}\ri{1};      
    for i = 2:level
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + RRi{i}\(ri{i} - APi{i}*ei{i});
    end
    eaux = [II*reshape(ei{level},dim*N,1); zeros(Ndof-size(IIt,2),1)]; % transfer back to the edge element space
    % Part2:  grad*(BP)^{-1}*grad^t
    % ri = cell(level,1);
    ri{level} = gradt*r(1:NE);          % transfer the residual to the null space    
    for i = level:-1:2
        ei{i} = Si{i}\ri{i};   
        ri{i-1} = Res{i}*(ri{i}-BPi{i}*ei{i});
    end
    ei{1} = BPi{1}\ri{1};      
    for i=2:level
        ei{i} = ei{i} + Pro{i-1}*ei{i-1};
        ei{i} = ei{i} + SSi{i}\(ri{i}-BPi{i}*ei{i});
    end
    eaux = eaux + [grad*ei{level}; zeros(Ndof-NE,1)]; % transfer back to the edge element space
    % combined correction in the finest grid
    Br = eh + eaux;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
