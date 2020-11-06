function [u,A,assembleTime,solverTime] = PoissonVEM(node,elem,pde,option)
%% POISSONVEM solve Poisson equation using virtual element method
%
%     -\Delta u = f,  in Omega
%             u = g_D on the boundary of Omega
%
% in a domain described by node and elem, with boundary edges Dirichlet.
% Each element could be different polygon shape.
%
% Input:
%   node, elem: coordiante and connvectivity;
%   pde.f: functional handle right side or data
%   pde.g_D: functional handle for Dirichelet condition
%
% Output:
%   u: solution on the current mesh
%   A: stiffness matrix
%
% Example
%
%   suarePoissonVEM
%
% See also Poisson, suarePoissonVEM
%
% The code is based on the following reference but optimized using
% vectorization to avoid for loops.
%
% Author: Min Wen and Long Chen
%
% Reference: Programming of Linear Virtual Element Methods. Long Chen and
% Min Wen. 2020. 
%
% Reference:  'The Hitchhiker's guide to the virtual element method'.
% by L.Beirao da Veiga, F.Brezzi, L.D.Marini, A.Russo.2013
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Preprocess
if ~exist('option','var'), option = []; end
N = size(node,1); % number of nodes
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isfield(option,'bdFlag'), option.bdFlag = []; end

%% Assemble the matrix equation
elemVertexNumber = cellfun('length',elem);% number of vertices per element
nnz = sum(elemVertexNumber.^2); % a upper bound on non-zeros
ii = zeros(nnz,1); % initialization
jj = zeros(nnz,1);
ss = zeros(nnz,1);
b = zeros(N,1);
edge = zeros(sum(elemVertexNumber),2);
index = 0;
edgeIdx = 1;
tic;
fv = pde.f(node); % right hand side evaluated at vertices
for nv = min(elemVertexNumber):max(elemVertexNumber)
    % find polygons with nv vertices
    idx = (elemVertexNumber == nv); % index of elements having nv vertices
    NT = sum(idx); % number of elements having nv vertices
    if NT == 0     % no element has nv vertices
        continue;
    end
    % vertex index and coordinates
    nvElem = cell2mat(elem(idx));
    x1 = reshape(node(nvElem,1),NT,nv);
    y1 = reshape(node(nvElem,2),NT,nv);
    x2 = circshift(x1,[0,-1]);
    y2 = circshift(y1,[0,-1]);
    % record edges
    nextIdx = edgeIdx + NT*nv;
    newEdgeIdx = edgeIdx:nextIdx-1;
    edge(newEdgeIdx,1) = nvElem(:); % get edge per element
    vertexShift = circshift(nvElem,[0,-1]);
    edge(newEdgeIdx,2) = vertexShift(:);
    edgeIdx = nextIdx;
    % compute geometry quantity: edge, normal, area, center
    bdIntegral = x1.*y2 - y1.*x2;
    area = sum(bdIntegral,2)/2; % the area per element
    h = repmat(sign(area).*sqrt(abs(area)),1,nv); % h is not the diameter
    cx = sum(reshape(node(nvElem(:),1),NT,nv),2)/nv;
    cy = sum(reshape(node(nvElem(:),2),NT,nv),2)/nv;
    normVecx = y2 - y1; % normal vector is a rotation of edge vector
    normVecy = x1 - x2;
    % matrix B, D, I - P
    Bx = (normVecx + circshift(normVecx,[0,1]))./(2*h); % average of normal vectors
    By = (normVecy + circshift(normVecy,[0,1]))./(2*h); % in adjaency edges
    Dx = (x1 - repmat(cx,1,nv))./h; %  m = (x - cx)/h
    Dy = (y1 - repmat(cy,1,nv))./h;
    IminusP = zeros(NT,nv,nv);
    for i = 1:nv
        for j = 1:nv
            IminusP(:,i,j) = - 1/nv - Dx(:,i).*Bx(:,j) - Dy(:,i).*By(:,j);
        end
        IminusP(:,i,i) = ones(NT,1) + IminusP(:,i,i);
    end
    % diffusion coefficient
    if ~isempty(pde.d) && isnumeric(pde.d)
        K = pde.d(idx);                                 % d is an array
    end
    if ~isempty(pde.d) && ~isnumeric(pde.d)       % d is a function
        K = pde.d([cx, cy]);
    end
    % assemble the matrix
    for i = 1:nv
        for j = 1:nv
            ii(index+1:index+NT) = nvElem(:,i);
            jj(index+1:index+NT) = nvElem(:,j);
            ss(index+1:index+NT) = K.*(Bx(:,i).*Bx(:,j) + By(:,i).*By(:,j) ...
                                 + dot(IminusP(:,:,i),IminusP(:,:,j),2));
            index = index + NT;
        end
    end
    % compute the right hand side
    patchArea = accumarray(nvElem(:),repmat(area/nv,nv,1),[N 1]); 
    b = b + fv.*patchArea;
end
A = sparse(ii,jj,ss,N,N);
assembleTime = toc;

%% Find boundary edges and nodes
totalEdge = sort(edge(:,1:2),2);
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
bdEdge  = [j(s==1), i(s==1)]; % find the boundary edge
isBdNode = false(N,1);
isBdNode(bdEdge) = true;

%% Impose Dirichlet boundary conditions
u = zeros(N,1);
u(isBdNode) = pde.g_D(node(isBdNode,:));
b = b - A*u;

%% TODO:modify Neuman boundary

%% Solve Au = b
isFreeNode = ~isBdNode; % all interior nodes are free
tic;
u(isFreeNode) = A(isFreeNode,isFreeNode)\b(isFreeNode);
solverTime = toc;