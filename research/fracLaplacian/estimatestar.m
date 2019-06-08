function eta = estimatestar(node,elem,u,pde,bdFlag)

N = size(node,1);

%% Assemble the equation for P2 element
option.solver = 'none';
[u2,Du,eqn] = PoissonP2(node,elem,pde,bdFlag,option);

%% Prolongate P1 solution to P2 space
NE = size(eqn.edge,1);
P1toP2 = sparse([(1:N)'; N+(1:NE)'; N+(1:NE)'], ...
                [(1:N)'; double(eqn.edge(:))],...
                [ones(N,1); 0.5*ones(2*NE,1)]',Nxdof,N);
u = P1toP2*u;

%% Local problems on stars
r = eqn.b - eqn.A*u;
isBdEdge = true(NE,1);
isBdEdge(eqn.freeDof(eqn.freeDof > N) - N) = false;
isBdNode = true(N,1);
isBdNode(eqn.freeDof(eqn.freeDof <= N)) = false;
intEdgeIdx = find(~isBdEdge);
e2v = sparse([intEdgeIdx,intEdgeIdx], double(eqn.edge(intEdgeIdx,:)), 1, NE, N);
estar = zeros(N,1);
for i = 1:N
    starEdges = find(e2v(:,i));
    if isBdNode(i)
        starDof = starEdges+N;
    else
        starDof = [i; starEdges+N];
    end
    localA = A(starDof(:),starDof(:));
    etastar = localA\r(starDof(:));
    estar(i) = etastar'*localA*etastar;
end

%% Change nodal wise to element wise indicator
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
estar = estar./valence;
eta = sqrt(estar(elem(:,1)) + estar(elem(:,2)) + estar(elem(:,3)));