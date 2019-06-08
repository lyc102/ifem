function [u,p,info] = asmgstokes(A,B,f,g,u,p,node,elem,bdFlag,ufreeDof,option)
%% ASMGSTOKES: auxiliary space multigird methods for Stokes equations
%
%  Created by Ming Wang, at April, 2012. Rewritten by Long Chen.

t = cputime;

%% Size of systems
Nu = length(u);
u = u(ufreeDof);
Np = length(p);
Ndof = length(f) + length(g);

%% Options
% Assign default values to unspecified parameters
if ~exist('option','var')
    option = []; 
end
inputoption = option; % record the input one
option = mgoptions(option,Ndof);    % parameters
% specific choices for mgstokesRT0
if ~isfield(inputoption,'smoothingstep')  % smoothing steps
    option.smoothingstep = 2;
end
if ~isfield(inputoption,'smootherSp')  % smoothing steps
    option.smootherSp = 'GS';
end
if ~isfield(inputoption,'smootherbarSp')  % smoothing steps
    option.smootherbarSp = 'SGS';
end
tol = option.tol; maxIt = option.solvermaxit; 
printlevel = option.printlevel;

%% Auxiliary Matrix for LSC-DGS smoother
Bt = B';
BBt = B*Bt;
BABt = B*A*Bt;
Su = tril(A);
Sp = tril(BBt);
Spt = triu(BBt);
DSp = diag(BBt);
auxMat = struct('Bt',Bt,'BBt',BBt,'BABt',BABt,'Su',Su,...
                   'Spt',Spt,'Sp',Sp,'DSp',DSp);
Ai_IU = []; Si_IU = []; SSi_IU = []; Res_IU = []; Pro_IU = [];
if strcmp(option.smootherbarSp,'VCYCLE')
%     option.smootherbarSp = 'VCYCLE';
    innerMGoption.solver = 'NO';
    innerMGoption.N0 = 10;
    [tempvar,tempvar,Ai_IU,Si_IU,SSi_IU,Res_IU,Pro_IU] = ...
                mg(auxMat.BBt,ones(Np,1),elem,innerMGoption); 
end

%% Assemble matrix of RT0-P0 elements
% data structure
elemunSort = elem;
[elem,bdFlag] = sortelem(elemunSort,bdFlag);
[elem2edge,edge] = dofedge(elem);
[Clambda,area,elemSign] = curlbasis(node,elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1); 

% matrices
vecMv = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
Mv = spdiags(vecMv,0,N,N);
Me = getmassmatvec(elem2edge,area,Clambda,'RT0');
invMt = spdiags(1./area,0,NT,NT);
B0 = -icdmat(double(elem2edge),elemSign*[1 -1 1]);
C = icdmat(double(edge),[-1 1]);
R = spdiags(1./vecMv,0,N,N)*C'*Me;
A0 = B0'*invMt*B0 + R'*Mv*R;

% Find Dirichlet boundary dof: fixedDof and pDof
isFixedEdge = false(NE,1);
if ~isempty(bdFlag) % if bdFlag is not empty
    isDirichlet(elem2edge(bdFlag(:)==1)) = true;
    isFixedEdge(isDirichlet) = true;% dof on D-edges
end
fixedEdge = find(isFixedEdge);
freeEdge = find(~isFixedEdge);
A0 = A0(freeEdge,freeEdge);
B0 = B0(:,freeEdge);

%% Transfer from TMAC to the current element
% velocity u
switch Nu
    case 2*(N+NE) % P2 or isoP2 element 
        pro_u = transferRT0toP2(node,elem,fixedEdge);
        pro_u = pro_u(ufreeDof,freeEdge);
        res_u = pro_u';
    case 2*NE  % CR element
        pro_u = transferRT0toCR(node,elem,fixedEdge);
        pro_u = pro_u(ufreeDof,freeEdge);
        res_u = pro_u';
    case 2*(N+NT) % Mini element
        pro_u = transferRT0toP1b(node,elem);
        pro_u = pro_u(ufreeDof,freeEdge);
        res_u = pro_u';
    case 2*NE + NT  % BDM1B element
        pro_u = speye(length(ufreeDof),length(freeEdge));
        res_u = pro_u';
end
% pressure p
switch Np
    case N  % P1 element
        [pro_p,res_p] = transferP0toP1(node,elem);
    case NT % P0 element
        pro_p = speye(NT);
        res_p = pro_p';
end

%% DCMG
% set up multilevel operators
optionRT = option;
optionRT.printlevel = 0;
optionRT.solver = 'NO';
zerou = zeros(NE,1);
zerop = zeros(NT,1);
[tempvaru,tempvarup,info,Ai,Bi,invMt,Pro_u,Pro_p] = ...
mgstokesRT0(A0,B0,zerou,zerop,zerou,zerop,node,elemunSort,freeEdge,optionRT);
% then set up parameters
optionRT.mu = 2;
optionRT.solver = 'WCYCLE';
optionRT.solvermaxit = 1;
optionRT.setupflag = false;
% optionRT.printlevel = 0;

nb = norm([f;g]);
err = zeros(maxIt,1);
err(1) = 1;
k = 1;
while (err(k) > tol) && (k <= maxIt)
    k = k+1;
    
    % pre-smoothing
    [u,p] = StokesLSCDGS(u,p,f,g,A,B,auxMat,elemunSort,option,...
                         Ai_IU,Si_IU,SSi_IU,Res_IU,Pro_IU);
    
    % form residual
    ru = f - A*u - Bt*p;
    rp = g - B*u;
    err(k) = norm([ru;rp])/nb;
    if printlevel >= 2    
        fprintf('#dof: %8.0u, iter: %2.0u, relres = %8.4e\n',...
                 Ndof, k, err(k));
    end
    % transfer to TMAC
    ru_RT0 = res_u*ru;
    rp_RT0 = res_p*rp;
    
    % solve with auxiliary space
    [eu_RT0,ep_RT0] = mgstokesRT0(A0,B0,ru_RT0,rp_RT0,zerou,zerop,...
                node,elemunSort,freeEdge,optionRT,Ai,Bi,invMt,Pro_u,Pro_p);

    % transfer back to the current element
    eu = pro_u*eu_RT0;
    ep = pro_p*ep_RT0;
    
    % correction on fine space
    u = u + eu;
    p = p + ep;
    
    %  post-smoothing
    [u,p] = StokesLSCDGS(u,p,f,g,A,B,auxMat,elemunSort,option,...
                         Ai_IU,Si_IU,SSi_IU,Res_IU,Pro_IU);
end
err = err(1:k);
itStep = k-1;

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
time = cputime - t;
if printlevel >= 1
    fprintf('#dof: %6.0u,  #nnz: %6.0u, MG %6s iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
             Ndof, nnz(A)+2*nnz(B), 'ASMG', itStep, err(end), time)    
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));