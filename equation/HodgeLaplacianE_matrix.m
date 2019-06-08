function [Abar,isFreeEdge,isFreeNode,Mv,G,C] = HodgeLaplacianE_matrix(node,elem,bdFlag)
%% HODEGELAPLACIANE_MATRIX get matrices of Hodge Laplacian of edge element
%
% [Abar,isFreeEdge,isFreeNode,Mv,G,C] = HodgeLaplacianE_matrix(node,elem,bdFlag)
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
[elem,bdFlag]    = sortelem(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dofedge(elem);
[Dlambda,area]   = gradbasis(node,elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Assemble matrix 
Mv = accumarray(elem(:),[area;area;area]/3,[N,1]);
invMv = spdiags(1./Mv,0,N,N);
Mv = spdiags(Mv,0,N,N); % lumped mass matrix is diagonal
Me = getmassmatvec(elem2edge,area,Dlambda,'ND0');
invMt = spdiags(1./area,0,NT,NT);
G  = Me*icdmat(double(edge),[-1 1]);  % gradient matrix
R  = icdmat(double(elem2edge),[1 -1 1]);
C  = R'*invMt*R;

%% Boundary conditions
if isempty(bdFlag)
    %Dirichlet boundary condition only
    bdFlag = setboundary(node,elem,'Dirichlet');
end
if ~isempty(bdFlag)
    % Find boundary edges and nodes
    isBdEdge = false(NE,1);
    isBdNode = false(N,1);
    isBdEdge(elem2edge(bdFlag(:) == 1)) = true;
    bdEdge = edge(isBdEdge,:);
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