function [w,u,eqn,info] = biharmonicP1(node,elem,bdFlag,pde,option)
%%
% Add by Jie Zhou. Clean up is needed.

%%
if nargin<5, option = []; end
N = size(node,1); NT = size(elem,1);
tic;

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
A = sparse(N,N);
B = sparse(N,N);
for i = 1:3
    for j = i:3   
        %%
        % $A_{ij}|_{\tau} = \int_{\tau}K \phi_i\cdot \phi_j dxdy$ 
          Aij=1/factorial(1+1+2)*2*area;       
          if(i==j)
          Aij=2*Aij;   % A is the mass matrix.
          end
          Bij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
            B = B + sparse(elem(:,i),elem(:,j),Bij,N,N);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Aij; Aij],N,N); 
            B = B + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Bij; Bij],N,N);                            
        end        
    end
end
clear K Aij Dphi Bij 
%% Assemble the right hand side  
b = zeros(N,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[N 1]);
end
clear pxy bt lambda area 


%%   Find Dirichlet boundary nodes and modify the stiffness matrix and right vector
%    Find Dirichlet boundary nodes: fixedNode
fixedNode = [];
if ~isempty(bdFlag) % bdFlag specifies different bd conditions
    allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
    Dirichlet = allEdge((bdFlag(:) == 1),:);
    isBdNode = false(N,1); 
    isBdNode(Dirichlet(:)) = true;
    fixedNode = find(isBdNode);
    freeNode = find(~isBdNode);
end
%

[b1,u] = getbdP1(b);

B(:,fixedNode)=[];
Nu=size(freeNode,1);             
b(fixedNode)=[];

%% Solve the system of linear equations
% Solver
switch solver
    case 'direct'
        t = cputime;
        bigA = [A, B; ...
                B', sparse(Nu,Nu)];
        bigF = [b1; -b];
        bigu = bigA\bigF;    
        w = bigu(1:N);
        u(freeNode) = bigu(N+1:end);
        residual = norm(bigF-bigA*bigu);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);                
end % end of four order.

%     soln = struct('u',u,'Du',Du);
    eqn = struct('A',A,'B',B','f',b1,'g',-b,'freeNode',freeNode,'Lap',A);
    info.assembleTime = assembleTime;

    function [b1,u] = getbdP1
    %% Boundary conditions for Poisson equation: P1  FEM.
    %  Modify the right hand side b. The Neumann boundary integral is added
    %  to b. For Dirichlet boundary ndoes, b(fixedDof) is the evaluation of
    %  pde.g_D.
    %   
    % The order of assigning Neumann and Dirichlet boundary condition is
    % important to get the right setting at the intersection nodes of Dirichlet
    % and Neumann boundary edges.

    u  = zeros(N,1);
%     b1 = zeros(N,1); 
    %% Part 1: Find boundary edges and modify the load b    
    % Neumann boundary condition
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;  
        end   

        [tempvar,Neumann] = findboundary(elem);
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),2);
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNu = pde.g_N(ppxy);
            for igN = 1:2
                ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNu;
            end
        end
        ge = ge.*repmat(el,1,2);
        b1 =  accumarray(Neumann(:), ge(:),[N,1]); 

    %% Part 2: Find Dirichlet boundary edges and compute the boundary value
    % Dirichlet boundary conditions
         u(fixedNode) = pde.g_D(node(fixedNode,:));
                   b1 = b1 - B*u;

    end % end of getbdP1
end