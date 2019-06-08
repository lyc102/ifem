function [u,w,AE,AI,isExteriorTElem,isExteriorSElem] = interfacefittedPoisson(node,telem,selem,pde,interfaceEdge,bdEdge,option)
%% INTERFACEFITTEDPOISSON Poisson equation: P1 linear element.
%
%   u = INTERFACEPOISSON(node,telem,selem,pde,E) produces the linear finite element
%   approximation of the interface Poisson equation
% 
%   input:
%       node, N*2 matrix, node(i,:) are the xy coordinates of i-th node;
%       telem, NTT*3 matrix, telem(j,:) are the three  global indices of the
%              vertices ofj-th triangle element;
%       selem,NTS*4 matrix, selem(j,:) are the four global indices of the
%             vertices of j-th quad element;
%
% 

[N,Dim] = size(node); 
Ndof = N;

%% Geometry structures of interface meshes

NTS = size(selem,1);
isExteriorSElem = false(NTS,1);
c = (node(selem(:,1),:) + node(selem(:,2),:) + node(selem(:,3),:)+node(selem(:,4),:))/4;
isExteriorSElem(pde.phi(c)>0) = true;

NTT = size(telem,1);
isExteriorTElem = false(NTT,1);
c = (node(telem(:,1),:) + node(telem(:,2),:) + node(telem(:,3),:))/3;
isExteriorTElem(pde.phi(c)>0) = true;

[sAE,sbE] = getstiffmatrixandrhsonquad(node,selem(isExteriorSElem,:));
[sAI,sbI] = getstiffmatrixandrhsonquad(node,selem(~isExteriorSElem,:));
[tAE,tbE] = getstiffmatrixandrhsontri(node,telem(isExteriorTElem,:));
[tAI,tbI] = getstiffmatrixandrhsontri(node,telem(~isExteriorTElem,:));


A = sAE + sAI + tAE + tAI;
b = sbE + sbI + tbE + tbI;

AI = sAI + tAI;
AE = sAE + tAE;

flux = getfluxconditiononinterface(node,interfaceEdge);
b = b - flux;



%% Extend w to the whole domain
isInterfaceNode = false(Ndof,1);
isInterfaceNode(interfaceEdge(:)) = true;
interfaceNode = find(isInterfaceNode);

isInNode = false(Ndof,1);
interiorSElem = selem(~isExteriorSElem,:);
interiorTElem = telem(~isExteriorTElem,:);
isInNode([interiorSElem(:);interiorTElem(:)]) = true;
inNode = find(isInNode & ~isInterfaceNode);

w = zeros(Ndof,1);
FI= zeros(Ndof,1);
w(interfaceNode) = pde.exactw(node(interfaceNode,:));
FI = FI - AI*w;
extensionoption.printlevel = 0;
w(inNode) = amg(AI(inNode, inNode),FI(inNode),extensionoption);

%% Dirichlet boundary condition
isBdNode = false(Ndof,1); 
isBdNode(bdEdge(:)) = true;
bdNode = find(isBdNode);

u = zeros(Ndof,1); 
u(bdNode) = pde.g_D(node(bdNode,:));
b = b - A*u+AI*w;
b(bdNode) = u(bdNode);

freeNode = find(~isBdNode);
bdidx = zeros(Ndof,1); 
bdidx(bdNode) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
AD = T*A*T + Tbd;

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
end

%%
eqn = struct('A',AD,'b',b,'freeNode',freeNode);

   %% Assemble stiffness matrix and rhs on quad mesh
    function [A,b] = getstiffmatrixandrhsonquad(node,elem)
        [NT,NV] = size(elem);
        % generate sparse pattern
        ii = zeros(10*NT,1); jj = zeros(10*NT,1);
        index = 0;
        for i = 1:4
            for j = i:4
                ii(index+1:index+NT) = double(elem(:,i));
                jj(index+1:index+NT) = double(elem(:,j));
                index = index + NT;
            end
        end
        % quadrature points
        if ~isfield(pde,'d'), pde.d = []; end
        if ~isfield(option,'dquadorder')
            option.dquadorder = 2;        % default order is exact for quadratic function
        end
        [pts, weight] = quadptsquad(option.dquadorder);
        nQuad = size(pts,1);
        % compute non-zeros
        sA = zeros(10*NT,nQuad);
        for p = 1:nQuad
            % Dphi at quadrature points
            [phi, Dphip, J] = quadbasis(node,elem,pts(p,:));
            index = 0;
            for i = 1:4
                for j = i:4
                    Aij = 0;
                    if isempty(pde.d) || isnumeric(pde.d)
                        Aij = Aij + weight(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
                    else
                        pxy = zeros(NT, Dim);
                        for ip = 1:Dim
                            xi = node(:,ip);
                            pxy(:,ip) = xi(elem)*phi;
                        end
                        Aij = Aij + weight(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.d(pxy);
                    end
                    if ~isempty(pde.d) && isnumeric(pde.d) % d is piecewise constant
                        Aij = pde.d.*Aij;
                    end
                    Aij = Aij.*J;
                    sA(index+1:index+NT,p) = Aij;
                    index = index + NT;
                end
            end
        end
        sA = sum(sA,2);
        % assemble the matrix
        diagIdx = (ii == jj);   upperIdx = ~diagIdx;
        A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
        AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
        A = A + AU + AU';
        clear Aij ii jj Dphip
        
        %% Assemble the right hand side
        b = zeros(Ndof,1);
        if ~isfield(option,'fquadorder')
            option.fquadorder = 3;   % default order
        end
        if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
            pde.f = [];
        end
        
        if ~isempty(pde.f)
            [pts, weight] = quadptsquad(option.fquadorder);
            nQuad = size(pts,1);
            bt = zeros(NT,NV);
            for p = 1:nQuad
                % quadrature points in the x-y coordinate
                [phi, ~, J] = quadbasis(node,elem, pts(p,:));
                pxy = zeros(NT, Dim);
                for i = 1:Dim
                    xi = node(:,i);
                    pxy(:,i) = xi(elem)*phi; % ? questionable
                end
                fp = pde.f(pxy);
                bt = bt + (weight(p)*fp.*J)*phi';
            end
            b = accumarray(elem(:),bt(:),[Ndof 1]);
        end
        
    end

    %% Assemble the stiffmatrix and right hand side on triangle mesh
    function [A,b] = getstiffmatrixandrhsontri(node,elem)
        NT = size(elem,1);
        % quadrature points
        center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
       % Diffusion coefficient
        if isfield(pde,'d') && ~isempty(pde.d)
            if isnumeric(pde.d)
                K = pde.d;                   % d is an array
            else                            % d is a function
                K = pde.d(center);
            end
        else
            K = [];
        end
        [Dphi,area] = gradbasis(node,elem);
        A = sparse(Ndof,Ndof);
        for i = 1:3
            for j = i:3
                Aij = (Dphi(:,1,i).*Dphi(:,1,j)+Dphi(:,2,i).*Dphi(:,2,j)).*area;
                if ~isempty(K)
                    Aij = K.*Aij;
                end
                if (j==i)
                    A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,N);
                else
                    A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],[Aij; Aij],Ndof,Ndof);
                end
            end
        end

        b = zeros(Ndof,1);
        if ~isfield(option,'fquadorder')
            option.fquadorder = 3;   % default order
        end
        if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
            pde.f = [];
        end
        if ~isempty(pde.f)
            [lambda,weight] = quadpts(option.fquadorder);
            nQuad = size(lambda,1);
            ft = zeros(NT,3);
            for p = 1:nQuad
                % quadrature points in the x-y coordinate
                pxy = lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:);
                % function values at quadrature points
                fp = pde.f(pxy);
                % evaluate fp outside.
                for j = 1:3
                    ft(:,j) = ft(:,j) + lambda(p,j)*weight(p)*fp;
                end
            end
            ft = ft.*[area,area,area];
            b = accumarray(elem(:),ft(:),[Ndof 1]);
        end  
    end

    %% Neumann boundary condition on interface edges
    function b = getfluxconditiononinterface(node,interfaceEdge)
        ve = node(interfaceEdge(:,1),:) - node(interfaceEdge(:,2),:);
        ve = [-ve(:,2), ve(:,1)];
        edgeLen = sqrt(sum(ve.^2,2));
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 5;   % default order
        end
        [lambda,weight] = quadpts1(option.gNquadorder);
        nQuad = length(weight);
        ge = zeros(size(interfaceEdge,1),2);
        for i = 1:nQuad
            pxy=lambda(i,1)*node(interfaceEdge(:,1),:)+lambda(i,2)*node(interfaceEdge(:,2),:);
            leftpt = pxy - ve/2;
            rightpt = pxy + ve/2;
            pxy = findintersectbisect(pde.phi,leftpt,rightpt);
            q = pde.exactq(pxy);
            ge(:,1)=ge(:,1)+weight(i)*lambda(i,1)*q;
            ge(:,2)=ge(:,2)+weight(i)*lambda(i,2)*q;
        end
        ge = ge.*[edgeLen,edgeLen];
        b = accumarray(interfaceEdge(:), ge(:),[Ndof,1]);
    end

end