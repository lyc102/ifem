function uh = CurlCurl2dNd0(node,elem,bdFlag,pde)
%% CURLCURL2DND0 curlcurl equation: lowest order edge element on triangular mesh
% rot(\mu^{-1} curl u) + \kappa u = f in \Omega
% u\times n = 0 on \Gamma_D
% \mu^{-1} curl u = g_N on \Gamma_N

% Solve u_h in Nd0 \cap {u \times n = 0 on \Gamma_D} = V_h
% (\mu^{-1} curl u_h, curl v_h) + (\kappa u_h, v_h)
% = (f, v_h) + <v_h \cdot \tau, g_N>_{\Gamma_N} for any v_h \in V_h

% basis function associated with an edge |e|
% \phi_e = \frac{|e|}{2|K|} (\lamda_{i+1} n_{i-1} - \lambda_{i+1} n_{i-1})
% implented based on an older iFEM circa 2010.
% 
% For 3D curl-curl problems
% See also Maxwell, cubeMaxwell
%
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
%%
NT = size(elem,1);
T = auxstructure(elem);
elem2edge = T.elem2edge;
edge = T.edge;
elem2edgeSign = ones(NT,3);
totalEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
idx = (totalEdge(:,1)>totalEdge(:,2));
elem2edgeSign(idx) = -1;

if ~isfield(pde,'mu'), pde.mu = 1; end
if ~isfield(pde,'kappa'), pde.kappa = 1; end
if ~isfield(pde,'g_D'), pde.g_D = []; end
if ~isfield(pde,'g_N'), pde.g_N = []; end

NE = size(edge,1);
uh = zeros(NE,1);
g_D = pde.g_D;
g_N = pde.g_N;
f_x = pde.f_x;
f_y = pde.f_y;
mu = pde.mu;
kappa = pde.kappa;


%% Boundary
idxD = (bdFlag(:) == 1);     % all Dirichlet edges in bdFlag
isFixedEdge = false(NE,1);
isFixedEdge(elem2edge(idxD)) = true;  % index of fixed boundary edges
freeEdge = ~isFixedEdge;

idxN = (bdFlag(:) == 2);     % all Neumann edges in bdFlag
Neumann = edge(elem2edge(idxN),:);

%% geometric quantities
%edge vector which follows elem2edge counterclockwisely
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

%length_ve(:,:,i) is the length of i^th edge
length_ve = sqrt(sum(ve.^2,2));
length_ve = reshape(length_ve,NT,3);

%outer normal vector of each edge which follows the index of elem2edge
ne = [ve(:,2,:), -ve(:,1,:)];


%%
%M is the mass matrix M_{ij} = (\kappa \phi_i, \phi_j)
%B is the stiffness matrix B_{ij} = (\mu^{-1} curl \phi_i, curl \phi_j)
c = [3 1 2 3 1];
%cyclic group of {1,2,3}: c(i)=i-1, c(i+2)=i+1

M = sparse(NE,NE);
B = sparse(NE,NE);
muinv = 1/mu;

for i = 1:3
    for j = 1:3
        sij = elem2edgeSign(:,i).*elem2edgeSign(:,j);
        
        if i==j
            nini = length_ve(:,i).^2;
            njnk = dot(ne(:,:,c(i)),ne(:,:,c(i+2)),2);
            Mij = kappa.*length_ve(:,i).*(nini - 3*njnk)./(24*area);
        else
            ninj = dot(ne(:,:,i),ne(:,:,j),2);
            ni2 = length_ve(:,i).^2;
            nj2 = length_ve(:,j).^2;
            Mij = -sij.*kappa.*length_ve(:,i).*length_ve(:,j).*(ni2 + nj2 + 3*ninj)./(24*area);
        end
        
        
        Bij = muinv.*sij.*length_ve(:,i).*length_ve(:,j)./area;
        
        M = M + sparse(elem2edge(:,i),elem2edge(:,j),Mij,NE,NE);
        B = B + sparse(elem2edge(:,i),elem2edge(:,j),Bij,NE,NE);
    end
end

A = M+B;

%% rhs (old hard-coded implementation)
Fx = quadelem2node(node,elem,f_x);
Fy = quadelem2node(node,elem,f_y);

bt1x = -0.5*elem2edgeSign(:,1).*length_ve(:,1)...
    .*(Fx(:,2).*ne(:,1,3) - Fx(:,3).*ne(:,1,2))./area;
bt2x = -0.5*elem2edgeSign(:,2).*length_ve(:,2)...
    .*(Fx(:,3).*ne(:,1,1) - Fx(:,1).*ne(:,1,3))./area;
bt3x = -0.5*elem2edgeSign(:,3).*length_ve(:,3)...
    .*(Fx(:,1).*ne(:,1,2) - Fx(:,2).*ne(:,1,1))./area;

bt1y = -0.5*elem2edgeSign(:,1).*length_ve(:,1)...
    .*(Fy(:,2).*ne(:,2,3) - Fy(:,3).*ne(:,2,2))./area;
bt2y = -0.5*elem2edgeSign(:,2).*length_ve(:,2)...
    .*(Fy(:,3).*ne(:,2,1) - Fy(:,1).*ne(:,2,3))./area;
bt3y = -0.5*elem2edgeSign(:,3).*length_ve(:,3)...
    .*(Fy(:,1).*ne(:,2,2) - Fy(:,2).*ne(:,2,1))./area;

bt1 = bt1x + bt1y;
bt2 = bt2x + bt2y;
bt3 = bt3x + bt3y;

F = accumarray(elem2edge(:),[bt1;bt2;bt3],[NE 1]);

%% Dirichlet Boundary
if ~isempty(pde.g_D)
    fixedEdge = edge(isFixedEdge,:);
    signFixedEdge = elem2edgeSign(idxD);
    Nve = node(fixedEdge(:,1),:) - node(fixedEdge(:,2),:);
    NveLength = sqrt(sum(Nve.^2,2)); 
    GD = quadedge(node,fixedEdge,g_D).*signFixedEdge./NveLength;
    uh(isFixedEdge) = GD;
end

%% Neumann Boundary
if (~isempty(Neumann) && ~isempty(pde.g_N))
    signNeumannEdge = elem2edgeSign(idxN);
    GN = quadedge(node,Neumann,g_N).*signNeumannEdge;
    F = F + accumarray(idxN(:),GN(:),[NE,1]); 
end

%direct solver
uh(freeEdge) = A(freeEdge,freeEdge)\F(freeEdge);
