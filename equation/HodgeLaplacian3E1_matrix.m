function [Abar,isFreeEdge,isFreeNode,Mv,G,C] = HodgeLaplacian3E1_matrix(node,elem,bdFlag)
%% HODEGELAPLACIAN3E1_MATRIX get matrices of Hodge Laplacian of the linear edge element
%
% [Abar,isFreeEdge,isFreeNode,Mv,G,C] = HodgeLaplacian3E1_matrix(node,elem,bdFlag)
% generates matrices of mixed finite element method of Hodge Laplacian
% using the linear order edge element space (Nd0 + grad(quadratic edge bubble)). 
%
% Output:
%  - Abar: Schur complement G*inv(Mv)*G' + R'*inv(Mt)*R;
%  - Mv: mass matrix
%  - G: gradient matrix
%  - C = R'*inv(Mt)*R
%  - isFreeEdge: logical index of interior edges (no boundary condition given)
%  - isFreeNode: logical index of interior nodes (no boundary condition given)
%  
%  Other key variables:
%  - isFreeDofS: logical index of free dofs for sigma, two copies of isFreeEdge 
%  - isFreeDofu: logical index of free dofs for u, freeDof for P2 elements
%
% The saddle point system is in the form 
%      |-M   G| |sigma|  = |f|
%      | G'  C| |u|      = |g|
%
% It is used in mgHodgeLapE to assemble matrices in the coase levels.
%
% See also mgHodgeLapE
% 
% Modified from HodgeLaplacian3E_matrix copying Long's routine for linear ND.
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
[elem,bdFlag]    = sortelem3(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dof3edge(elem);
[Dlambda,volume] = gradbasis3(node,elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NdofU = N+NE; NdofS = 2*NE;
%% Assemble matrix 
% Mass matrices
elem2dofu = dof3P2(elem);
Mv = getmassmat3(node,elem2dofu,volume);
MvLump = diag(Mv);
invMv = spdiags(1./MvLump,0,NdofU,NdofU);
Mv = spdiags(MvLump,0,NdofU,NdofU);
Me = getmassmatvec3(elem2edge,volume,Dlambda,'ND1');
grad = icdmat(double(edge),[-1,1]); % P1 to ND0
% e2v matrix: P1 to P2 lambda_i -> lambda_i(2lambda_i - 1)
ii = [(1:NE)';(1:NE)'];
jj = double(edge(:));
ss = -2*[ones(NE,1),ones(NE,1)];
e2v = sparse(ii,jj,ss,NE,N);
bigGrad = [grad, sparse(NE,NE); ...
            e2v, 4*speye(NE,NE)];
G = Me*bigGrad;

% R'MR: curlcurl operator
%locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
curlPhi(:,:,6) = 2*mycross(Dlambda(:,:,3),Dlambda(:,:,4),2);
curlPhi(:,:,1) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,2),2);
curlPhi(:,:,2) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,3),2);
curlPhi(:,:,3) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,4),2);
curlPhi(:,:,4) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,3),2);
curlPhi(:,:,5) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,4),2);
ii = zeros(21*NT,1); jj = zeros(21*NT,1); sC = zeros(21*NT,1); 
index = 0;
for i = 1:6
    for j = i:6
        % local to global index map
        % curl-curl matrix
        Cij = dot(curlPhi(:,:,i),curlPhi(:,:,j),2).*volume;
        ii(index+1:index+NT) = double(elem2edge(:,i)); 
        jj(index+1:index+NT) = double(elem2edge(:,j));
        sC(index+1:index+NT) = Cij;
        index = index + NT;
    end
end
clear curlPhi % clear large size data
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
C = sparse(ii(diagIdx),jj(diagIdx),sC(diagIdx),NdofS,NdofS);
CU = sparse(ii(upperIdx),jj(upperIdx),sC(upperIdx),NdofS,NdofS);
C = C + CU + CU';

%% Boundary conditions
if isempty(bdFlag)
    %Dirichlet boundary condition only
    bdFlag = setboundary(node,elem,'Dirichlet');
end
if ~isempty(bdFlag)
    % Find boundary edges and nodes
    isBdEdge = false(NE,1);
    isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
    isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
    isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
    isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
    bdEdge = edge(isBdEdge,:);
    isBdNode = false(N,1);
    isBdNode(bdEdge(:)) = true;
    isFreeEdge = ~isBdEdge;
    isFreeNode = ~isBdNode;
    isFreeDofS = [isFreeEdge; isFreeEdge];
    isFreeDofu = [isFreeNode; isFreeEdge];
end

%% Restrict the matrix to free variables
G  = G(isFreeDofS,isFreeDofu);
C  = C(isFreeDofS,isFreeDofS);
Mv = Mv(isFreeDofu,isFreeDofu);
% Schur complement
Abar  = G*invMv(isFreeDofu,isFreeDofu)*G' + C;