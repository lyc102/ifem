function [u,p,info] = mgDarcy(M,B,f,g,elem,option,varargin)
%% MGDARCY multigrid solvers for Darcy system discretized by RT0-P0 element
%
% [u,p,info] = mgDarcy(M,B,f,g,elem) solves saddle point problem
% discretized from RT0-P0 mixed FEM for Darcy equation.
%
%      |M  B'| |u|  = |f|
%      |B  0 | |p|  = |g|
%
%  A V-cycle multigrid using overlapping Schwarz smoother is implemented.
%  In the first step, mgdivDarcy is called to find initial u satisfying Bu
%  = g. Then at each vertex patch, a local problem with prescribed boundary
%  flux is solved. Detailed description and convergence analysis can be
%  found in the reference.
%
%  It is around 10 times slower than tripremixPoisson since a large for
%  loop is not efficient in MATLAB. In operation count, this is superior.
%  See the complexity analysis in page 18 of the reference.
%
% Reference: L. Chen. Multigrid Methods for Constrained Minimization
% Problems and Application to Saddle Point Problems
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

time = cputime;
%% Size of systems
Nu = length(f);                    
Np = length(g);
N = max(elem(:));                  % number of nodes
% NT = size(elem,1);               % number of elements

%% Options
% Assign default values to unspecified parameters
if ~exist('option','var'), 
    option = []; 
end
option = mgoptions(option,Nu);    % parameters
u0 = option.x0; 
N0 = option.N0; 
tol = option.tol;
maxIt = option.solvermaxit; 
mu = option.smoothingstep; 
coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; 
freeEdge = option.freeEdge;

%% Hierarchical Structure of Mesh
HB = zeros(N,3);
level = 20;
NL(level+1) = N; % now NL(1:level) = 0;
elemi = cell(level,1);
elemi{level} = elem;
freeEdgei = cell(level,1);
freeEdgei{level} = freeEdge;
Pro_u = cell(level,1);
Pro_p = cell(level,1);
for k = level: -1 : 2
    switch option.refType 
        case 'red'
            [elemi{k-1},newHB] = uniformcoarsenred(elemi{k}); % coasen red refinement
            if ~isempty(newHB)
                [Pro_u{k-1}, freeEdgei{k-1}] = transferedgered(elemi{k-1},elemi{k},freeEdgei{k}); % transfer operator of u
                Pro_p{k-1} = repmat(speye(size(elemi{k-1},1)),4,1);   % transfer operator of p
            end
        case 'bisect'
            % merge two coarsen of bisection grids s.t. the ratio is 1/4.
            [tempelem,newHB1,tree] = uniformcoarsen(elemi{k}); % coarse bisection
            if ~isempty(newHB1)            
                % first coarsen
                [tempPro_u,tempFreeEdge] = transferedgecoarsen(tempelem,elemi{k},tree,freeEdgei{k});
                tempPro_p = transferelem(tempelem,elemi{k},tree);
                % second coarsen
                [elemi{k-1},newHB2,tree] = uniformcoarsen(tempelem); % coarse bisection
                [Pro_u{k-1}, freeEdgei{k-1}] = transferedgecoarsen(elemi{k-1},tempelem,tree,tempFreeEdge);
                Pro_u{k-1} = tempPro_u*Pro_u{k-1};
                Pro_p{k-1} = transferelem(elemi{k-1},tempelem,tree);
                Pro_p{k-1} = tempPro_p*Pro_p{k-1};
                newHB = [newHB1; newHB2];
            end
    end    
    if (isempty(newHB)) || (size(elemi{k},1)< 2*N0) 
        % no nodes are removed or it reaches the coarsest level
        NL = NL(k:end);       
        break; 
    end
    NL(k) = NL(k+1) - size(newHB,1); % update NL(k)
    HB(NL(k)+1:NL(k+1),1:3) = newHB(:,1:3);
end
level = length(NL)-1;    % actual level
elemi = elemi(end-level+1: end);
Pro_u = Pro_u(end-level+1: end);
Pro_p = Pro_p(end-level+1: end);
freeEdgei = freeEdgei(end-level+1: end);
% generate edge, elem2edge and freeNode etc
edgei = cell(level,1);
elem2edgei = cell(level,1);
for k = 1:level
    [elem2edgei{k},edgei{k}] = dofedge(elemi{k});
end

%% Matrices in each level
oldf = f;
Mi = cell(level,1);
Bi = cell(level,1);
if size(M,1) > length(freeEdge) % truncate to free edge only
    Mi{level} = M(freeEdge,freeEdge);    
    Bi{level} = B(:,freeEdge);
    f  = f(freeEdge);
    u0 = u0(freeEdge);
else
    Mi{level} = M;    
    Bi{level} = B;    
end
Res_u = cell(level,1);
Res_p = cell(level,1);
for k = level:-1:2
    Res_u{k} = transpose(Pro_u{k-1});
    Res_p{k} = transpose(Pro_p{k-1});
    Mi{k-1} = Res_u{k}*Mi{k}*Pro_u{k-1};           
    Bi{k-1} = Res_p{k}*Bi{k}*Pro_u{k-1};
end

% %% Exact solver: for debug only
% C = sparse(size(B,1),size(B,1));
% A = [M B';B C];
% F = [f; zeros(size(B,1),1)];
% tempu = A\F;
% exactSigma = tempu(1:Ndof);

%% MG cycles
% initial set up
k = 1; 
f0 = f - Mi{level}*u0;
g0 = g - Bi{level}*u0;
% find u satisfy Bu = g
u = u0 + mgdivDarcy(f0,g0);
% temp = abs(g-Bi{level}*u);
% idx = temp<1e-14;
% findelem(node,elem,idx,'noindex','FaceColor','c');  
if printlevel >= 1
    fprintf('Multigrid Vcycle Iteration \n')
end
err = zeros(maxIt,1);
err(1) = 1;
r0 = f - Mi{level}*u;
r = r0;
while (max(err(k,:)) > tol) && (k <= maxIt)
    k = k + 1;
    % Step 2: Compute Br by one Vcylce MG
    Br = vcycle(r);
    % Step 3: Correct the solution
    u = u + Br;
    % Step 1: Form residual r
    r = r - Mi{level}*Br;
    err(k) = sqrt(abs(Br'*r/(u'*r0))); % approximate relative error in energy norm
    if printlevel >= 2
        fprintf('#dof: %8.0u,  #nnz: %8.0u, MG Vcycle iter: %2.0u, err = %8.4e\n',...
                 Nu+N, nnz(M)+nnz(B), k-1, err(k));
    end            
end

err = err(1:k,:);
itStep = k-1;

%% Find pressure
Ap = B*B';
ubar = zeros(Nu,1);
ubar(freeEdge) = u;
b = B*(oldf-M*ubar);
p = zeros(Np,1);
if option.isPureNeumannBC
    freep = 1:Np-1;
else
    freep = 1:Np;
end
p(freep) = Ap(freep,freep)\b(freep);

%% Output
time = cputime - time;
if k > maxIt
    flag = 1;
else
    flag = 0;
end
if printlevel >= 2
    fprintf('#dof: %8.0u, level: %2.0u,   coarse grid %2.0u, #nnz: %8.0u\n',...
              Nu+N, level, size(Mi{1},1), nnz(Mi{1}))
end
if printlevel >= 1
    fprintf('#dof: %8.0u,  #nnz: %8.0u, smoothing: %2.0u, iter: %2.0u,   err = %8.4e,   time = %4.2g s\n',...
                 Nu+Np, nnz(M)+2*nnz(B), mu, itStep, err(end), time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions vcycle
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
        ei{i} = zeros(size(Mi{i},1),1);
        ei{i} = SchwarzsmootherDarcy(Mi{i},Bi{i},ri{i},zeros(size(Bi{i},1),1),...
                ei{i},elemi{i},edgei{i},elem2edgei{i},freeEdgei{i},mu);
        ri{i-1} = Res_u{i}*(ri{i} - Mi{i}*ei{i});
    end
    if strcmp(coarsegridsolver,'direct')
        Ncoarse = size(Bi{1},1);
        C1 = sparse(Ncoarse,Ncoarse);
        A1 = [Mi{1} Bi{1}'; Bi{1} C1];
        F1 = [ri{1}; zeros(Ncoarse,1)];
        if option.isPureNeumannBC
            bigu = A1(1:end-1,1:end-1)\F1(1:end-1); % pure Neumann
        else
            bigu = A1\F1;
        end
        ei{1} = bigu(1:size(Mi{1},1));
    end
    for i = 2:J
        ei{i} = ei{i} + Pro_u{i-1}*ei{i-1};
        ei{i} = SchwarzsmootherDarcy(Mi{i},Bi{i},ri{i},zeros(size(Bi{i},1),1),...
                ei{i},elemi{i},edgei{i},elem2edgei{i},freeEdgei{i},-mu);
    end
    Br = ei{J};
    end
%% div Darcy
% Use one V-cycle with post-smoothing only to find u s.t. div u = g

    function u = mgdivDarcy(rf,rg)
    J = level;
    rfi = cell(J,1);
    rgi = cell(J,1);
    eui = cell(J,1);
    rfi{J} = rf;
    rgi{J} = rg;
    % restrict the residual to the coarse level
    for i = J:-1:2
        rfi{i-1} = Res_u{i}*rfi{i};
        rgi{i-1} = Res_p{i}*rgi{i};
    end
    % exact solve in the coarest mesh
    Ncoarse = size(Bi{1},1);
    C1 = sparse(Ncoarse,Ncoarse);
    A1 = [Mi{1} Bi{1}'; Bi{1} C1];
    F1 = [rfi{1}; rgi{1}];
    if option.isPureNeumannBC
        bigu = A1(1:end-1,1:end-1)\F1(1:end-1); % pure Neumann
    else
        bigu = A1\F1;
    end
    eui{1} = bigu(1:size(Mi{1},1));
    % prolongate to the fine level and post-smoothing in each element
    for i = 2:J
        eui{i} = Pro_u{i-1}*eui{i-1};
        eui{i} = SchwarzsmootherelemDarcy(Mi{i},Bi{i},rfi{i},rgi{i},...
                eui{i},elemi{i},edgei{i},elem2edgei{i},freeEdgei{i},1);
    end
    u = eui{J};
    end
end