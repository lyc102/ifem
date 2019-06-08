function [w,u,eqn,info] = biharmonicP2(node,elem,bdFlag,pde,option)
%%
%   [w,u] = biharmonicP1(node,elem,bdFlag,pde) produces the mixed cubic finite element
%   approximation of the biharmonic equation, where w = laplace u
%   See also biharmonicP1, biharmonicP2, biharmonicP3.
%
%   Created by Jie Zhou. Clean up is needed
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin<5, option = []; end
tic;
%% Construct Data Structure
[elem2dof,edge,bdDof] = dofP2(elem);
N = size(node,1);  NT = size(elem,1);  Ndof = N+size(edge,1);

%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
% Since Dphi_i*Dphi_j is quadratic, numerical quadrature rule is used here
if ~isfield(option,'quadorder')
    option.quadorder = 4;   % default order
end
[lambda, weight] = quadpts(option.quadorder);
nQuad = size(lambda,1);
ii = zeros(21*NT,1); jj = zeros(21*NT,1); sA = zeros(21*NT,1);sB = zeros(21*NT,1);
index = 0;
for i = 1:6
    for j = i:6
        Bij = 0;
        Aij = 0;        
        for p = 1:nQuad
                 Bij = Bij + weight(p)*dot(Dphi(p,i),Dphi(p,j),2);
                 Aij = Aij + weight(p)*dot(phi(p,i),phi(p,j),2);
        end
        Bij = Bij.*area;
        Aij = Aij.*area;

        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));
        sB(index+1:index+NT) = Bij;
        sA(index+1:index+NT) = Aij;        
        index = index + NT;
    end
end
clear Aij Bij
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
B = sparse(ii(diagIdx),jj(diagIdx),sB(diagIdx),Ndof,Ndof);
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
% A = spdiags(accumarray(ii(diagIdx),sA(diagIdx),[Ndof 1]),0,Ndof,Ndof);
BU = sparse(ii(upperIdx),jj(upperIdx),sB(upperIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
B = B + BU + BU';
A = A + AU + AU';

%% boundary condition
%%
    fixedDof = [];
    isFixedDof = false(Ndof,1); 
    if ~isempty(bdFlag)     
        elem2edge = elem2dof(:,4:6)-N;
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        isFixedDof(edge(isDirichlet,:)) = true;
        isFixedDof(N + find(isDirichlet')) = true;
        fixedDof = find(isFixedDof);
        freeDof = find(~isFixedDof);    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction Dphi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function s = Dphi(p,i) % gradient of basis phi
    switch i
        case 1
            s = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
        case 2
            s = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
        case 3
            s = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
        case 4
            s = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        case 5
            s = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
        case 6
            s = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction Dphi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function s = phi(p,i) % gradient of basis phi
    switch i
        case 1
            s = (2*lambda(p,1)-1).*lambda(p,1);           
        case 2
            s = (2*lambda(p,2)-1).*lambda(p,2);            
        case 3
            s = (2*lambda(p,3)-1).*lambda(p,3);             
        case 4
            s = 4*lambda(p,3).*lambda(p,2); 
        case 5
            s = 4*lambda(p,1).*lambda(p,3); 
        case 6
            s = 4*lambda(p,1).*lambda(p,2); 
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assemble right hand side by high order quadrature rule
% To reduce the effect of the error introduced by the numerical quadrature,
% the load term is computed using the 3rd order qudrature rule.
b = zeros(Ndof,1);
% u=zeros(Ndof,1);
w = zeros(Ndof,1);

if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
%     phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
%     phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
%     phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
%     phi(:,4) = 4*lambda(:,2).*lambda(:,3);
%     phi(:,5) = 4*lambda(:,3).*lambda(:,1);
%     phi(:,6) = 4*lambda(:,1).*lambda(:,2);
    bt = zeros(NT,6);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        if isfield(pde,'f') && isnumeric(pde.f)
            fp = pde.f;        % piecewise constant       
        else
            fp = pde.f(pxy);   % function handle
        end
        for j = 1:6
            bt(:,j) = bt(:,j) + w(p)*phi(p,j)*fp;
        end
    end
    bt = bt.*repmat(area,1,6);
    b = accumarray(elem2dof(:),bt(:),[Ndof 1]); 
end





[b1,u] = getbdP2(b);
B(:,fixedDof)=[];
Nu=size(freeDof,1);

b(fixedDof)=[];

switch solver
    case 'direct'
        t = cputime;
        bigA = [A, B; ...
                B', sparse(Nu,Nu)];
        bigF = [b1; -b];
        bigu = bigA\bigF;    
        w = bigu(1:Ndof);
        u(freeDof)=bigu(Ndof+1:Ndof+Nu);
        residual = norm(bigF-bigA*bigu);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);                
end % end of four order.

%     soln = struct('u',u,'Du',Du);
    eqn = struct('A',A,'B',B','f',b1,'g',-b,'freeNode',freeNode,'Lap',A);
    info.assembleTime = assembleTime;


 function [b1,u] = getbdP2
    %% Boundary conditions for Poisson equation: P2 quadratic FEM.
    %
    % The set up of boundary condition consists of two parts: 
    %

    %
    %  Modify the right hand side b. The Neumann boundary integral is added
    %  to b. For Dirichlet boundary ndoes, b(fixedDof) is the evaluation of
    % pde.g_D.
    %
    u = zeros(Ndof,1);
    
    %% Part 1: Find boundary edges and modify the load b    
    % Neumann boundary condition
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 6;  
        end   
                idxN = (bdFlag(:) == 1);      % all Neumann edges in bdFlag        
        Neumannidx = elem2edge(idxN ); % index of Neumann and Robin edges
        % since boundary integral is also needed for Robin edges
        Neumann   = edge(Neumannidx,:);
        b1 = zeros(Ndof,1);
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        % quadratic bases (1---3---2)
        bdphi = zeros(nQuadgN,3);        
        bdphi(:,1) = (2*lambdagN(:,1)-1).*lambdagN(:,1);
        bdphi(:,2) = (2*lambdagN(:,2)-1).*lambdagN(:,2);
        bdphi(:,3) = 4*lambdagN(:,1).*lambdagN(:,2);
        % length of edge
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        ge = zeros(size(Neumann,1),3);
        for pp = 1:nQuadgN
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNu = pde.g_N(ppxy);
            ge(:,1) = ge(:,1) + weightgN(pp)*gNu*bdphi(pp,1);
            ge(:,2) = ge(:,2) + weightgN(pp)*gNu*bdphi(pp,2);
            ge(:,3) = ge(:,3) + weightgN(pp)*gNu*bdphi(pp,3); % interior bubble
        end
        % update RHS
        ge = ge.*repmat(el,1,3);        
        b1(1:N) = accumarray(Neumann(:), [ge(:,1); ge(:,2)],[N,1]);
        b1(N+Neumannidx) = b1(N+Neumannidx) + ge(:,3);

    %% Part 2: Find Dirichlet boundary edges and compute the boundary value
    % Dirichlet boundary conditions
   
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        % interpolation
                idx = (fixedDof > N);                         % index of edge nodes
        u(fixedDof(~idx)) = pde.g_D(node(fixedDof(~idx),:)); % bd value at vertex dofs
        bdEdgeIdx = fixedDof(idx) - N;
        bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
        u(fixedDof(idx)) = pde.g_D(bdEdgeMid);
        b1 = b1 - B*u;
  

    end % end of getbdP2

end                 % end of function biharmonicP2