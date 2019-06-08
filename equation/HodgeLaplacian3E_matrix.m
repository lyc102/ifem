function [Abar,isFreeEdge,isFreeNode,Mv,G,C] = HodgeLaplacian3E_matrix(node,elem,bdFlag)
%% HODEGELAPLACIAN3E_MATRIX get matrices of Hodge Laplacian of edge element
%
% [Abar,isFreeEdge,isFreeNode,Mv,G,C] = HodgeLaplacian3E_matrix(node,elem,bdFlag)
% generates matrices of mixed finite element method of Hodge Laplacian
% using the lowest order edge element space. 
%
% Output:
%  - Abar: Schur complement G*inv(Mv)*G' + R'*inv(Mt)*R;
%  - Mv: mass matrix
%  - G: gradient matrix
%  - C = R'*inv(Mt)*R
%  - isFreeEdge: logic index of free edges (no boundary condition given)
%  - isFreeNode: logic index of free nodes (no boundary condition given)
%
% The saddle point system is in the form 
%      |-M   G| |sigma|  = |f|
%      | G'  C| |u|      = |g|
%
% It is used in mgHodgeLapE to assemble matrices in the coase levels.
%
% See also mgHodgeLapE
% 
% Created by Jie Zhou and revised by Long Chen.
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
[elem,bdFlag]    = sortelem3(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dof3edge(elem);
[Dlambda,volume] = gradbasis3(node,elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Assemble matrix 
% Mass matrices
Mv = accumarray(elem(:),[volume;volume;volume;volume]/4,[N,1]);
invMv = spdiags(1./Mv,0,N,N);
Mv = spdiags(Mv,0,N,N);
Me = getmassmatvec3(elem2edge,volume,Dlambda,'ND0');
G = Me*icdmat(double(edge),[-1 1]);  % gradient matrix
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
C = sparse(ii(diagIdx),jj(diagIdx),sC(diagIdx),NE,NE);
CU = sparse(ii(upperIdx),jj(upperIdx),sC(upperIdx),NE,NE);
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
end

%% Restrict the matrix to free variables
G  = G(isFreeEdge,isFreeNode);
C  = C(isFreeEdge,isFreeEdge);
Mv = Mv(isFreeNode,isFreeNode);
% Schur complement
Abar  = G*invMv(isFreeNode,isFreeNode)*G' + C;