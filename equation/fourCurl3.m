function [w,u,eqn,info] = fourCurl3(node,elem,bdFlag,pde,option)
%% FOURCURL3: the fourth order curl problem in 3D
%
% [w,u,eqn,info] = fourCurl3(node,elem,bdFlag,pde,option)
% uses the lowest order Nedelec elements to approximate the velocity u and
% the stream function w. 
%
% We solve the following equations:
%
%  -w + curl curl u = 0 
%  curl curl w + u  = f
%
% with Dirichlet boundary condition
%
% u\times n = (curl u )\times n  = 0  on \partial \Omega
%
% Please check fourCurl3doc for details.
%
% Lin Zhong, May, 2013.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


if ~exist('option','var'), option = []; end

%% Data structure 
elemunSort = elem;
[elem,bdFlag] = sortelem3(elemunSort,bdFlag);
[elem2edge,edge] = dof3edge(elem);
[elem2face,face] = dof3face(elem);

NE = size(edge,1);
NF = size(face,1);
NT = size(elem,1);

t = cputime;

%% Assemble Matrix
[Dlambda,volume] = gradbasis3(node,elem);

% Mass matrix in H(curl)
Me = getmassmatvec3(elem2edge,volume,Dlambda,'ND0');

% Mass matrix in H(div)
Mf = getmassmatvec3(elem2face,volume,Dlambda,'RT0');

% incidence matrix of vertex-edge
face2edge = zeros(NF,3,'int32');
face2edge(elem2face(:,1),:) = elem2edge(:,[4 5 6]);
face2edge(elem2face(:,2),:) = elem2edge(:,[2 3 6]);
face2edge(elem2face(:,3),:) = elem2edge(:,[1 3 5]);
face2edge(elem2face(:,4),:) = elem2edge(:,[1 2 4]);
C = icdmat(double(face2edge),[1,-1,1]);

% Whole matrix
B = C'*Mf*C;
A = [-1*Me B; B Me];

%% Righthand side
locEdge = [1 2; 1 3; 1 4; 2 3;2 4; 3 4];
f = zeros(NE,1);

if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end

if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end

if isfield(pde,'f') && ~isempty(pde.f)
    [lambda,w] = quadpts3(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,6);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ... 
            + lambda(p,3)*node(elem(:,3),:)...
            + lambda(p,4)*node(elem(:,4),:);
        fp = pde.f(pxyz);
        for k = 1:6
            i = locEdge(k,1); j = locEdge(k,2);
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            phi_k = lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i);
            rhs = dot(phi_k,fp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end    
    bt = bt.*repmat(volume,1,6);
    f = accumarray(elem2edge(:),bt(:),[NE 1]);
end
clear pxyz fp bt rhs phi_k

bigf = [zeros(NE,1);f];
assembleTime = cputime - t;

%% Boundary condition
isBdEdge = [];
if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N)
    % Dirichlet boundary condition only
    bdFlag = setboundary3(node,elem,'Dirichlet');
end

if ~isempty(bdFlag)
    % Find boundary edges and nodes
    isBdEdge = false(NE,1);
    isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
    isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
    isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
    isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
end
freeEdge = ~isBdEdge;

if any(freeEdge)
    idx = [true(NE,1); freeEdge];
    A_bd= A(idx,idx);
    bigf_bd = bigf(idx);
end

%% Solve by direct solver
t = cputime;
bigu = A_bd\bigf_bd;
w = bigu(1:NE);
u = zeros(NE,1);
u(freeEdge) = bigu(NE+1:end);
info.solverTime = cputime -t;
display(info.solverTime);

%% Output information
eqn = struct('elem2edge',elem2edge,'elem2face',elem2face,'B',B,'edge',edge,'Me',Me,'f',f);
info.assembleTime = assembleTime;

end