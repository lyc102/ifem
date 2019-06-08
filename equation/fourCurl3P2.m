function [omega,u,eqn,info] = fourCurl3P2(node,elem,bdFlag,pde,option)

%% fourCurl3P2: the fourth order curl problem in 3D
%
% [w,u,eqn,info] = fourCurl3P2(node,elem,bdFlag,pde,option)
% uses the first P2 elements to approximate the velocity u and
% the stream function w. 
%
% We solve the following equations:
%  -w + curl curl u = 0 
%  curl curl w + u  = f
% u\times n = (curl u )\times n  = 0  on \partial \Omega
%
% please check fourCurl3doc for details.
%
% Lin Zhong, June 2013.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Set up optional input arguments
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

%% Sort elem to asend ordering
elemunSort = elem;
[elem,bdFlag] = sortelem3(elemunSort,bdFlag);

%% Construct Data Structure
[Dlambda,volume] = gradbasis3(node,elem);
volume = abs(volume);
%---------------------------------------------%
% locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % phi
               1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % psi
               3 2 4; 3 1 4; 2 1 4; 2 1 3; ...
               4 2 3; 4 1 3; 4 1 2; 3 1 2]; % chi
           
%---------------------------------------------%
[elem2edge,edge] = dof3edge(elem);
[elem2face,face] = dof3face(elem);
%---------------------------------------------%
NE = size(edge,1);
NF = size(face,1);
NT = size(elem,1);
Ndof = 2*(NE+NF);
%---------------------------------------------%
elem2dof = [elem2edge elem2edge+NE elem2face+2*NE elem2face+2*NE+NF];
face2edge = zeros(NF,3,'int32');
face2edge(elem2face(:,1),:) = elem2edge(:,[4 5 6]);
face2edge(elem2face(:,2),:) = elem2edge(:,[2 3 6]);
face2edge(elem2face(:,3),:) = elem2edge(:,[1 3 5]);
face2edge(elem2face(:,4),:) = elem2edge(:,[1 2 4]);

tic;

%% Assemble Matrix
DiDjcross = zeros(NT,3,4,4);
for i = 1:4
    for j = i+1:4        
        DiDjcross(:,:,i,j) = mycross(Dlambda(:,:,i),Dlambda(:,:,j),2);
        DiDjcross(:,:,j,i) = -DiDjcross(:,:,i,j);
    end
end
DiDj = zeros(NT,4,4);
for i = 1:4
    for j = i:4        
        DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end

ii = zeros(210*NT,1); jj = zeros(210*NT,1); 
index = 0;
for i = 1:20
    for j = i:20
        ii(index+1:index+NT) = double(elem2dof(:,i));
        jj(index+1:index+NT) = double(elem2dof(:,j));
        i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
        j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
        Aij = zeros(NT,1);  Mij = zeros(NT,1); 
        %% curl-curl matrix
        if (i<=6) && (j<=6)
            Aij = 4*dot(DiDjcross(:,:,i1,i2),DiDjcross(:,:,j1,j2),2);
        end
        if (i<=6) && (j>12)
            Aij = dot(DiDjcross(:,:,i1,i2),...
            DiDjcross(:,:,j1,j3)-DiDjcross(:,:,j1,j2)+2*DiDjcross(:,:,j2,j3),2)/2;
        end
        if (i>12) && (j>12)
        % curl chi_i =  Dlambda_{i1}mycross phi_{i2,i3} + lambda_{i1}curl phi_{i2,i3}
        % curl chi_j =  Dlambda_{j1}mycross phi_{j2,j3} + lambda_{j1}curl phi_{j2,j3}
        % (Dlambda_{i1}mycross phi_{i2,i3}) dot (Dlambda_{j1}mycross phi_{j2,j3})
            temp11 = ((1+(i2==j2))*dot(DiDjcross(:,:,i1,i3),DiDjcross(:,:,j1,j3),2) ...
                    - (1+(i2==j3))*dot(DiDjcross(:,:,i1,i3),DiDjcross(:,:,j1,j2),2) ...
                    - (1+(i3==j2))*dot(DiDjcross(:,:,i1,i2),DiDjcross(:,:,j1,j3),2) ...
                    + (1+(i3==j3))*dot(DiDjcross(:,:,i1,i2),DiDjcross(:,:,j1,j2),2))/20;
            % lambda_{i1}curl phi_{i2,i3} dot lambda_{j1}curl phi_{j2,j3}
            temp22 = (1+(i1==j1))/5*dot(DiDjcross(:,:,i2,i3),DiDjcross(:,:,j2,j3),2);
            % Dlambda_{i1}mycross phi_{i2,i3} dot lambda_{j1}curl phi_{j2,j3}
            temp12 = dot( (1+(j1==i2))*DiDjcross(:,:,i1,i3) ...
                        - (1+(j1==i3))*DiDjcross(:,:,i1,i2), ...
                                       DiDjcross(:,:,j2,j3),2)/10;                         
            % Dlambda_{j1}mycross phi_{j2,j3} dot lambda_{i1}curl phi_{i2,i3}
            temp21 = dot( (1+(i1==j2))*DiDjcross(:,:,j1,j3) ...
                        - (1+(i1==j3))*DiDjcross(:,:,j1,j2), ...
                                       DiDjcross(:,:,i2,i3),2)/10;   
            Aij = temp11 + temp22 + temp12 + temp21;
        end
        if (i<=6) && (j<=6)
            % block 1: (phi_i,phi_j)
            Mij = 1/20*((1+(i1==j1))*DiDj(:,i2,j2) ...
                      - (1+(i1==j2))*DiDj(:,i2,j1) ...
                      - (1+(i2==j1))*DiDj(:,i1,j2) ...
                      + (1+(i2==j2))*DiDj(:,i1,j1));
        end
        if (i<=6) && (7<=j) && (j<=12)
            % block 2: (psi_j,phi_i)
            Mij = 1/20*( (1+(j1==i1))*DiDj(:,j2,i2) ...
                       - (1+(j1==i2))*DiDj(:,j2,i1) ...
                       + (1+(j2==i1))*DiDj(:,j1,i2) ...
                       - (1+(j2==i2))*DiDj(:,j1,i1));

        end
        if (7<=i) && (i<=12) && (7<=j) && (j<=12)
            % block 3: (psi_j,psi_i)
            Mij = 1/20*((1+(i1==j1))*DiDj(:,i2,j2) ...
                      + (1+(i1==j2))*DiDj(:,i2,j1) ...
                      + (1+(i2==j1))*DiDj(:,i1,j2) ...
                      + (1+(i2==j2))*DiDj(:,i1,j1));
        end
        if (i<=6) && (j>12)
            % block 4: (chi_j,phi_i)
            Mij = intlambda([j1,i1,j2],3)*DiDj(:,i2,j3) ...
                 -intlambda([j1,i1,j3],3)*DiDj(:,i2,j2) ...
                 -intlambda([j1,i2,j2],3)*DiDj(:,i1,j3) ...
                 +intlambda([j1,i2,j3],3)*DiDj(:,i1,j2);            
        end
        if (7<=i) && (i<=12) && (j>12)
            % block 5: (chi_j,psi_i)
            Mij = intlambda([j1,i1,j2],3)*DiDj(:,i2,j3) ...
                 -intlambda([j1,i1,j3],3)*DiDj(:,i2,j2) ...
                 +intlambda([j1,i2,j2],3)*DiDj(:,i1,j3) ...
                 -intlambda([j1,i2,j3],3)*DiDj(:,i1,j2);                        
        end
        if (i>12) && (j>12)
            % block 6: (chi_j,chi_i)
            Mij = intlambda([i1,j1,i2,j2],3)*DiDj(:,i3,j3) ...
                 -intlambda([i1,j1,i2,j3],3)*DiDj(:,i3,j2) ...
                 -intlambda([i1,j1,i3,j2],3)*DiDj(:,i2,j3) ...
                 +intlambda([i1,j1,i3,j3],3)*DiDj(:,i2,j2);                        
        end
        Aij = Aij.*volume;
        Mij = Mij.*volume;
        sA(index+1:index+NT) = Aij;
        sM(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
clear curlBasis_i curlBasis_j basis_i basis_j
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
M = M + MU + MU';
clear AU MU

% Whole matrix
bigA = [-M A; A M];

%% Righthand side
f = zeros(Ndof,1);
if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   % default order is 4
end
if isfield(pde,'f') && ~isempty(pde.f)
    [lambda,w] = quadpts3(option.fquadorder); % quadrature order is 4
    nQuad = size(lambda,1);
    bt = zeros(NT,20);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:) ... 
             + lambda(p,4)*node(elem(:,4),:);
        fp = pde.f(pxyz);
        for k = 1:20
            k1 = locBasesIdx(k,1); 
            k2 = locBasesIdx(k,2); 
            k3 = locBasesIdx(k,3);
            % evaluate basis at quadrature points
            if k<=6
            % phi_k = lambda_{k1}Dlambda_{k2} - lambda_{k2}Dlambda_{k1};
                basis_k = (lambda(p,k1)*Dlambda(:,:,k2) ...
                          -lambda(p,k2)*Dlambda(:,:,k1));
            elseif k<=12
            % psi_k = lambda_{k1}Dlambda_{k2} + lambda_{k2}Dlambda_{k1};
                basis_k = (lambda(p,k1)*Dlambda(:,:,k2) ...
                          +lambda(p,k2)*Dlambda(:,:,k1));
            else
            % chi_k = lambda_{k1}phi_{k2,k3};    
                basis_k = lambda(p,k1)*(lambda(p,k2)*Dlambda(:,:,k3) ...
                                       -lambda(p,k3)*Dlambda(:,:,k2));
            end            
            rhs = dot(basis_k,fp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(volume,1,20);
    f = accumarray(elem2dof(:),bt(:),[Ndof 1]);
end
clear pxy fp bt rhs basis_k

bigf = [zeros(Ndof,1); f];
assembleTime = toc;

%% Boundary condition
% Find Dirichlet boundary dof: fixedDof
if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N)
    % Dirichlet boundary condition only
    bdFlag = setboundary3(node,elem,'Dirichlet');
end
isBdDof = false(Ndof,1);

if ~isempty(bdFlag)
    %% Dirichlet boundary condition on edge dofs
    % Find boundary faces, edges and nodes
    isBdFace = false(NF,1);
    isBdFace(elem2face(bdFlag(:,1) == 1,1)) = true;
    isBdFace(elem2face(bdFlag(:,2) == 1,2)) = true;
    isBdFace(elem2face(bdFlag(:,3) == 1,3)) = true;
    isBdFace(elem2face(bdFlag(:,4) == 1,4)) = true;
%     boundaryFace = face(isBdFace,:);
    bdFace2edge = face2edge(isBdFace,:);
    isBdEdge = false(NE,1);
    isBdEdge(bdFace2edge(:)) = true;
    edgeBdDof = [find(isBdEdge); NE + find(isBdEdge)];
%     bdEdge = edge(isBdEdge,:);
%     isBdNode(bdEdge) = true;
%     bdNode = find(isBdNode);
    faceBdDof = 2*NE + [find(isBdFace); NF+find(isBdFace)];
%     edgeIdxMap = zeros(NE,1);
%     edgeIdxMap(isBdEdge) = 1:size(bdEdge,1);
%     bdFace2edge = edgeIdxMap(bdFace2edge);
    isBdDof(edgeBdDof) = true;
    isBdDof(faceBdDof) = true;
end

freeDof = ~isBdDof;

if any(freeDof)
    idx = [true(Ndof,1); freeDof];
    bigA_bd= bigA(idx,idx);
    bigf_bd = bigf(idx);
end

%% Solve
t = cputime;
bigu = bigA_bd\bigf_bd;
omega = bigu(1:Ndof);
u = zeros(Ndof,1);
u(freeDof) = bigu(Ndof+1:end);
info.solverTime = cputime -t;
display(info.solverTime);

%% Output information
eqn = struct('elem2edge',elem2edge,'elem2face',elem2face,'face2edge',face2edge,'A',A,'edge',edge,'face',face,'M',M,'f',f);
info.assembleTime = assembleTime;
end