function [u,A,assembleTime,solverTime] = PoissonVEM(node,elem,pde)
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
%     node = [0,0; 0,0.5; 0,1; 0.5,0;1,0;0.5,0.5;0.5,1;1,0.5,1,1];
%     elem = [1,4,6,2; 6,8,9,7;4,5,8,6;2,6,7,3];
%     pde.f = inline('ones(size(p),1)','p');
%     pde.exactu = inline('(-p(:,1).^2-p(:,2).^2)/4','p');
%     [u,A] = PoissonVEM(node,elem,pde);
%     uI = pde.exactu(node);
%     error = sqrt((u-uI)'*A*(u-uI));
%
% See also Poisson, PoissonS
%
% The code is based on the following reference but optimized using
% vectorization to avoid for loops.
%
% Reference:
%
% Reference:  'The Hitchhiker's guide to the virtual element method'.
% by L.Beirao da Veiga, F.Brezzi, L.D.Marini, A.Russo.2013
%
% Author: Min Wen and Long Chen
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Assemble the matrix equation
N = size(node,1); % The number of nodes
elemVertexNumber = cellfun('length',elem);% the number of vertices per element
nnz = sum(elemVertexNumber.^2);
ii = zeros(nnz,1); %initialization
jj = zeros(nnz,1);
ss = zeros(nnz,1);
b = zeros(N,1);
edge = zeros(sum(elemVertexNumber),2);
index = 0;
edgeIdx = 1;
tic;
fv = pde.f(node);
for nv = min(elemVertexNumber):max(elemVertexNumber)
    % find polygons with Nv vertices
    idx = (elemVertexNumber == nv); % index of elements with Nv vertices
    NT = sum(idx); % the number of elements having the same number of vertices
    if NT == 0     % no element has Nv vertices
        continue;
    end
    % vertex index and coordinates: vertices are counter-clockwise
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
    h = repmat(sign(area).*sqrt(abs(area)),1,nv); % h = sqrt(area) not the diameter
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
    % assemble the matrix
    for i = 1:nv
        for j = 1:nv
            ii(index+1:index+NT) = nvElem(:,i);
            jj(index+1:index+NT) = nvElem(:,j);
            ss(index+1:index+NT) = Bx(:,i).*Bx(:,j) + By(:,i).*By(:,j) ...
                                 + dot(IminusP(:,:,i),IminusP(:,:,j),2);
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

%% Solve Au = b
isFreeNode = ~isBdNode; % all interior nodes are free
tic;
u(isFreeNode) = A(isFreeNode,isFreeNode)\b(isFreeNode);
solverTime = toc;