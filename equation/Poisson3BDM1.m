function [u,sigma,A] = Poisson3BDM1(node,elem,pde,bdFace,solver,varargin)
%% POISSON3BDM1 Poisson equation: linear BDM element in 3-D.
%
% Example
%    example3BDM1
%
% Created by Ming Wang at Dec 28, 2010. Debuged by Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
[elem2dof,dofSign,face] = dof3BDM1(elem);
NT = size(elem,1); NF = size(face,1);
locFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2]; % normal vector is oriented outside the tetrahedral.
locBasesIdx = [locFace(:,[1 2 3]);locFace(:,[2 3 1]);locFace(:,[3 1 2])];
[Dlambda,volume] = gradbasis3(node,elem);

%% Compute element-wise basis.
% on each face, the divergence of the three bases are the same.
divPhi(:,4) = dot(Dlambda(:,:,1),mycross(Dlambda(:,:,3),Dlambda(:,:,2),2),2);
divPhi(:,1) = dot(Dlambda(:,:,2),mycross(Dlambda(:,:,3),Dlambda(:,:,4),2),2);
divPhi(:,2) = dot(Dlambda(:,:,1),mycross(Dlambda(:,:,4),Dlambda(:,:,3),2),2);
divPhi(:,3) = dot(Dlambda(:,:,1),mycross(Dlambda(:,:,2),Dlambda(:,:,4),2),2);
divPhi = repmat(divPhi,1,3);

%% Assemble matrix 
NdofSigma = 3*NF;
Ndofu = NT;
% Part: A. Mass matrix for RT0 element
A = sparse(NdofSigma,NdofSigma);
for i = 1:12
    for j = i:12 
% Local basis
		% local to global index map and its sign
		ii = double(elem2dof(:,i));
		jj = double(elem2dof(:,j));
        signij = double(dofSign(:,i).*dofSign(:,j));
        i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
        j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
        % computation of mass matrix --- (phi_i, phi_j) 
		Aij = 1/20*volume.*((1+(i1==j1))* ... 
              dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                  mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2));
        Aij = Aij.*signij; 
        if (j==i)
            A = A + sparse(ii,jj,Aij,NdofSigma,NdofSigma);
        else
            A = A + sparse([ii;jj],[jj;ii],[Aij; Aij],NdofSigma,NdofSigma);        
        end        
% Local basis
    end
end
clear Aij
% Part: B. divergence operator
B = sparse(Ndofu,NdofSigma);
for i = 1:12  
	Bi = divPhi(:,i).*volume.*double(dofSign(:,i));
	B = B + sparse((1:Ndofu),double(elem2dof(:,i)),Bi,Ndofu,NdofSigma);
end
C = sparse(Ndofu,Ndofu);
bigA = [A B';B C];
%% Assemble right hand side.% if ~isempty(bdEdge)
%     idx = (bdEdge(:,1) == 1);
%     %sigma(elem2dof(:,1)) = a integeral to compute the value;
% end

F = zeros(Ndofu,1);
if ~isempty(pde.f)
	[lambda,weight] = quadpts3(2);
	nQuad = size(lambda,1);
	for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ...
			+ lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxy);
		F = F - fp*weight(p).*volume;
	end
end
bigF((NdofSigma+1):(NdofSigma+Ndofu),1) = F;
clear lambda weight nQuad fp area

%% Boundary Conditions 
% if ~isempty(bdEdge)
%     idx = (bdEdge(:,1) == 1);
%     %sigma(elem2dof(:,1)) = a integeral to compute the value;
% end
% *Dirichlet boundary condition*
% T = auxstructure(elem);
% edge2elem = T.edge2elem;
% isInEdge = edge2elem(:,1)~=edge2elem(:,2); % index of interior edge.
% isBdEdge= ~isInEdge; 
% DbdEdge = find(isBdEdge); % Dirichlet boundary edges.
% % Compute the integration of g_D on the boundary.
% G = zeros(length(T.bdEdge),1);
% G(DbdEdge) = double(dofSign(edge2elem(DbdEdge,1)+NT*(edge2el what is the difference ? Why in this way?
%% Local basisem(DbdEdge,3)-1))) ...
%              .*g_D(1/2*(node(edge(DbdEdge,1),:)+node(edge(DbdEdge,2),:)));
% F1 = zeros(NdofSigma+Ndofu,1);
% F1(DbdEdge) = G(DbdEdge);
% F1(NdofSigma+1:Ndofu+NdofSigma) = F;

%% Solver
%set solver type
if (nargin<=4)
    if (NT < 1e3)    % Direct solver for small size systems
        solver = 'direct';
    else            % Multigrid-type solver for large size systems
        solver = 'UzawaPCG';
    end
end
if strcmp(solver,'notsolve');
    sigma=[]; u =[];
elseif strcmp(solver,'direct');
   % bigu(freeNode) = bigA(freeNode,freeNode)\bigF(freeNode);
bigu= bigA\bigF;
    sigma = bigu(1:3*NF);
    u = bigu(3*NF+1:end);
else strcmp(solver,'UzawaPCG');
    [sigma,u] = uzawapcg(bigA(1:3*NF,1:3*NF),bigA(3*NF+1:end,1:3*NF),bigF(1:3*NF),bigF(3*NF+1:end));
end
%% TODO: 1. Handle the Dirichlet boundary conditon. 2. Write M-lint
