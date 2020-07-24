function [soln,eqn,info] = StokesBDM1B(node,elem,bdFlag,pde,option)
%% STOKESRT0 Stokes equations: the lowest order BDM element in 2D.
%
%  [soln,eqn,info] = StokesBDM1B(node,elem,bdFlag,pde,option)
%  uses the lowest order of BDM1 element with bubble function
%  (BDM1B) to approximate velocity u and piecewise
%  constant element (P0)  to approximate pressure p, repectively.
%
%  We solve the following equation:
%       - grad div u + curl rot u + grad p  = f   in \Omega   
%                                 - div u   = 0   in \Omega   
%                                       u   = g   on \Gamma   
%
% Based on a version by Ming Wang. Revised by Lin Zhong. Discussed with Jie
% Zhou and Long Chen. Further clean up by Long Chen. Add MG by Long.

if ~exist('option','var'), option = []; end

%% Mesh and data structure
elemunSort = elem;
[elem,bdFlag] = sortelem(elemunSort,bdFlag);
[elem2edge,edge] = dofedge(elem);
[Clambda,area,elemSign] = curlbasis(node,elem);
elem2dof = [elem size(node,1) + elem2edge]; % P2 elements
edge = double(edge); 
elem2edge = double(elem2edge);

%% Assembling Matrix 
N = size(node,1); NT = size(elem,1); NE = size(edge,1); 
Nw = N+NE+NT; Nu = 2*NE+NT; Np = NT; Ndof = Nu + Np;

t = cputime;
% Mvb: Mass matrix with lumping for vertex: P2 element + bubble 
Mvb = getmassmatrixP2(elem2dof,area,'NBB');
invMvb = spdiags(1./diag(Mvb),0,Nw,Nw);
Mv = getmassmatrixP2(elem2dof,area,'NB');

% Meb: Mass matrix for edge: BDM1B element
Meb =  getmassmatvec(elem2edge,area,Clambda,'BDM1B');

% invMt: the inverse of Mass matrix for P0 element
invMt = spdiags(1./area,0,NT,NT);

% B: -divergence operator
RT0B = -icdmat(double(elem2edge),elemSign*[1 -1 1]);% -div for RT0
B = [RT0B sparse(NT,NE+NT)];

% C: curl operator
% curl for RT0
RT0C = icdmat(double(edge),[-1 1]); 

% curl in HB (hierat = cputime bases)
HC = blkdiag(RT0C,4*speye(NE,NE),speye(NT,NT)); 

% -------------   The transform matrix between bases  -------------
% NB                                             HB
% 1) (2\lambda_i -1)\lambda_i +1/9w_b            1) \lambda_i
% 2) 4\lambda_i \lambda_j -4/9w_b         <===>  2) 4\lambda_i\lambda_j
% 3) 27\lambda_1\lambda_2\lambda_3               3) 27\lambda_1\lambda_2\lambda_3  
%
% NB = HB*T, where T is the tranfer matrix
%   T = | 1      0    0 |
%       |1/2     1    0 |
%       |1/9    4/9   1 |
N2NE = sparse([edge(:,1);edge(:,2)],repmat(1:NE,1,2),1,N,NE);
N2NT = sparse([elem(:,1);elem(:,2);elem(:,3)], repmat(1:NT,1,3),1,N,NT);
NE2NT = sparse([elem2edge(:,1);elem2edge(:,2);elem2edge(:,3)],repmat(1:NT,1,3),1,NE,NT);
T = [speye(N,N) sparse(N,NE) sparse(N,NT);
    -1/2*N2NE'  speye(NE,NE)   sparse(NE,NT);
    1/9*N2NT' -4/9*NE2NT' speye(NT,NT)];
C = HC*T;
  
% R: weak rot operator
R = invMvb*C'*Meb;

% Vector Laplacian
A = B'*invMt*B + R'*Mvb*R;

%% Assemble right hand side
localEdge = [2,3; 1 3; 1,2];
fu = zeros(Nu,1);% the right hand side of u
g = zeros(Np,1); % the right hand side of p
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   % default order is 3
end
if ~isempty(pde.f) 
    elem2dofu = [elem2edge NE+elem2edge 2*NE+(1:NT)'];
    % quadrature points in the barycentric coordinate
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    ft = zeros(NT,7);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fup = pde.f(pxy);   
        for j = 1:3
            j1 = localEdge(j,1); j2 = localEdge(j,2);
            phi_j = lambda(p,j1)*Clambda(:,:,j2)-lambda(p,j2)*Clambda(:,:,j1);
            ft(:,j)   = ft(:,j)  + w(p)*dot(fup,phi_j,2);
            psi_j = lambda(p,j1)*Clambda(:,:,j2)+lambda(p,j2)*Clambda(:,:,j1);
            ft(:,3+j) = ft(:,3+j)+ w(p)*dot(fup,psi_j,2);
        end
        chi = 27*(lambda(p,1)*lambda(p,2)*Clambda(:,:,3) + ...
                  lambda(p,2)*lambda(p,3)*Clambda(:,:,1) + ...
                  lambda(p,3)*lambda(p,1)*Clambda(:,:,2));
        ft(:,7) = ft(:,7) + w(p)*dot(fup,chi,2);
    end
    ft = ft.*repmat(area,1,7);
    fu = accumarray(elem2dofu(:),ft(:),[Nu 1]);
end
clear pxy fup ft phi_j psi_j chi

%% Boundary condition and graddiv part of A
[u,p,ufreeDof,pDof,utbd] = getbdStokesBDM1B;
assembleTime = cputime - t;

%% Solve the system of linear equations
% set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 1e5  % Direct solver for small size systems
        solver = 'direct';
    else             % Multigrid-type  solver for large size systems
        solver = 'mg';
    end
else
    solver = option.solver;
end
% solve the system
t = cputime;
% get submatrices of ufreeDof
A0 = A(ufreeDof,ufreeDof);
B0 = B(:,ufreeDof);
f0 = fu(ufreeDof);
g0 = g;
switch solver
    case 'direct'
    bigA = [A0, B0'; ...
            B0, sparse(Np,Np)];
    bigF = [f0; g0];
    bigu = [u; p];
    bigFreeDof = [ufreeDof; Nu+pDof];
    bigu(bigFreeDof) = bigA(1:end-1,1:end-1)\bigF(1:end-1);
    u = bigu(1:Nu);
    p = bigu(Nu+1:end);
    residual = norm(bigF - bigA*[u(ufreeDof); p]);
    info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);        
    case 'mg'
        option.solver  = 'WCYCLE';
        [u(ufreeDof),p,info] = mgstokesBDM(A0,B0,f0,g0,u,p,node,elemunSort,ufreeDof,option);
    case 'asmg'
        [u(ufreeDof),p,info] = asmgstokes(A0,B0,f0,g0,u,p,node,elemunSort,bdFlag,ufreeDof,option);                              
end

%% Post-process
if length(pDof)~= Np
    c = sum(p.*area)/sum(area);
    p = p - c;
end
w = R*u + invMvb*utbd;
% drop the bubble part in the output
% u = u(1:2*NE);
w = w(1:N+NE);

%% Output
soln = struct('u',u,'p',p,'w',w);
eqn = struct('A',A0,'B',B0,'Me',Meb,'Mv',Mv,'f',f0,'g',g0,...
             'edge',edge,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;

%%======================================================================
% subfunction getbdStokesBDM1B
%%======================================================================
    function [u,p,ufreeDof,pDof,utbd] = getbdStokesBDM1B
        %% Boundary condtion of Stokes equations: BDM1B-P0 elements
        
        % Initial set up
        % Nw = N+NE+NT; Nu = 2*NE+NT; Np = NT; Ndof = Nw+Nu+Np;
        utbd = zeros(Nw,1); 
        u = zeros(Nu,1);
        p = zeros(Np,1);
        ufreeDof = (1:Nu)';
        pDof = (1:Np-1)';
        
        if ~exist('bdFlag','var'), bdFlag = []; end
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'g_R'), pde.g_R = []; end        
        if isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
            bdFlag = [];
        end

        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedEdge = false(NE,1);
        if ~isempty(bdFlag) % if bdFlag is not empty
            % find out the Dirichlet boundary dof
            isDirichlet(elem2edge(bdFlag(:)==1)) = true;
            Dirichlet = edge(isDirichlet,:);
            isFixedEdge(isDirichlet) = true;% dof on D-edges
            fixedEdge = find(isFixedEdge);
            fixedDof = [fixedEdge; fixedEdge+NE];
            ufreeDof = setdiff((1:Nu)',fixedDof);
            % construct the edgeSign
            edgeSign = ones(NE,1);
            idx = (bdFlag(:,1) ~= 0 ) & (elemSign == -1) ; % the first edge is on boundary
            edgeSign(elem2edge(idx,1)) = -1;
            idx = (bdFlag(:,2) ~= 0 ) & (elemSign ==  1) ; % the second edge is on boundary
            edgeSign(elem2edge(idx,2)) = -1;
            idx = (bdFlag(:,3) ~= 0 ) & (elemSign == -1) ; % the third edge is on boundary
            edgeSign(elem2edge(idx,3)) = -1;     
        end
        
        % Compute the boundary integral
        if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D ==0))
            % 1. Normal component of u is imposed strongly
            if (isnumeric(pde.g_D) && length(pde.g_D) == Nu)
                u(fixedDof) = pde.g_D(fixedDof);
            else
                u(fixedDof) = faceinterpolate(pde.g_D, node,edge(fixedEdge,:),'BDM1');
            end
            % 2. Tangential component of u is imposed weakly
            % 4th order for the line integral
            [lambdagD,wgD] = quadpts1(4);
            nQuadgD = size(lambdagD,1);
            % quadrat = cputime bases 1--3--2
            bdphi = zeros(nQuadgD,3);
            bdphi(:,1) = (2*lambdagD(:,1)-1).*lambdagD(:,1);
            bdphi(:,2) = (2*lambdagD(:,2)-1).*lambdagD(:,2);
            bdphi(:,3) = 4*lambdagD(:,1).*lambdagD(:,2);
            ve = node(Dirichlet(:,2),:) - node(Dirichlet(:,1),:);
            ge = zeros(size(Dirichlet,1),3);
            int_left = zeros(size(Dirichlet,1),2);
            int_right = zeros(size(Dirichlet,1),2);
            int_mid = zeros(size(Dirichlet,1),2);
            for pp = 1:nQuadgD
                ppxy = lambdagD(pp,1)*node(Dirichlet(:,1),:) ...
                     + lambdagD(pp,2)*node(Dirichlet(:,2),:);
                gDp = pde.g_D(ppxy);
                int_left = int_left + wgD(pp)*gDp*bdphi(pp,1);
                int_right = int_right + wgD(pp)*gDp*bdphi(pp,2);
                int_mid = int_mid + wgD(pp)*gDp*bdphi(pp,3); % interior bubble
            end 
            ge(:,1) = dot(int_left,ve,2).*edgeSign(fixedEdge);
            ge(:,2) = dot(int_right,ve,2).*edgeSign(fixedEdge);
            ge(:,3) = dot(int_mid,ve,2).*edgeSign(fixedEdge);       
            utbd(1:N) = accumarray(Dirichlet(:), [ge(:,1); ge(:,2)],[N,1]);
            utbd(N+fixedEdge) = ge(:,3);
        end
        
        % Modify the right hand side
        fu = fu- A*u - Meb*(C*(invMvb*utbd));
         g = g - B*u;
         g = g - mean(g);
    end
end