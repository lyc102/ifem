function [soln,eqn,info] = Poisson3CR(node,elem,bdFlag,pde,option,varargin)
%% POISSON3CR Poisson equation: Crouzeix-Raviart element in 3-D.
%
%
% Example
%
%   Poisson3CRfemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
NT = size(elem,1);

%% Diffusion coefficient
t = cputime;  % record assembling time
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                   % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d)
    [lambda,weight] = quadpts3(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,1);
    for p = 1:nQuad
	   pxyz = lambda(p,1)*node(elem(:,1),:) ...
		 	+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:) ...
            + lambda(p,4)*node(elem(:,4),:);
        K = K + weight(p)*pde.d(pxyz);           % d is a function   
   end
end

%% Construct data structure 
[elem2face,face] = dof3face(elem); 
[Dphi,volume] = gradbasis3(node,elem);
NF = size(face,1);
Ndof = NF;

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof); 
for i = 1:4
    for j = i:4
        % local to global index map
		ii = double(elem2face(:,i));
		jj = double(elem2face(:,j));
        % local stiffness matrix.
        Aij = 9*dot(Dphi(:,:,i),Dphi(:,:,j),2).*volume;
        if ~isempty(pde.d)
            Aij = K.*Aij;
        end   
        if(j==i)
         A = A + sparse(ii,jj,Aij,Ndof,Ndof);
        else
         A = A + sparse([ii;jj],[jj;ii],[Aij; Aij],Ndof,Ndof);
        end
    end
end 
clear K Aij

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts3(option.fquadorder);
    phi = 1-3*lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,4);
    for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:)...
             + lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxyz);
        for i = 1:4
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(volume,1,4);
    b = accumarray(elem2face(:),bt(:),[Ndof 1]);
end
clear pxyz bt

%% Get boundary condition 
[AD,b,u,freeFace,isPureNeumann] = getbd3CR(b);

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeFace), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
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
    u(freeFace) = AD(freeFace,freeFace)\b(freeFace);         
    residual = norm(b - AD*u);
    info = struct('solverTime',toc,'itStep',0,'error',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'error',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        if nargin>=6
            HB = varargin{1};
        else
            HB = [];
        end
        [u,info] = mg(AD,b,elem,option,face,HB);
    case 'amg'
        option.solver = 'CG';
        [u(freeFace),info] = amg(AD(freeFace,freeFace),b(freeFace),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    patchVolume = accumarray(elem2face(:),repmat(volume/4,4,1), [NF 1]); 
    uc = sum(u.*patchVolume)/sum(volume);
    u = u - uc;   % int u = 0
end

%% Compute Du
dudx = u(elem2face(:,1)).*Dphi(:,1,1)+u(elem2face(:,2)).*Dphi(:,1,2) ...
     + u(elem2face(:,3)).*Dphi(:,1,3)+u(elem2face(:,4)).*Dphi(:,1,4);
dudy = u(elem2face(:,1)).*Dphi(:,2,1)+u(elem2face(:,2)).*Dphi(:,2,2) ...
     + u(elem2face(:,3)).*Dphi(:,2,3)+u(elem2face(:,4)).*Dphi(:,2,4);
dudz = u(elem2face(:,1)).*Dphi(:,3,1)+u(elem2face(:,2)).*Dphi(:,3,2) ...
     + u(elem2face(:,3)).*Dphi(:,3,3)+u(elem2face(:,4)).*Dphi(:,3,4);
Du = -3*[dudx, dudy, dudz]; 

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'face',face,'freeFace',freeFace,'Lap',A);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd3CR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeFace,isPureNeumann] = getbd3CR(b)
%% GETBD3CR Boundary conditions for Poisson equation: Crouzeix-Raviart element in 3D.

    u = zeros(Ndof,1);
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end
 
    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    Robin = [];
    idxR = (bdFlag(:) == 3);      % index of Robin faces in bdFlag
    if any(idxR)    
        isRobin = false(Ndof,1);
        isRobin(elem2face(idxR)) = true;
        Robin = face(isRobin,:);  % Robin faces  
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        v12 = node(Robin(:,2),:)-node(Robin(:,1),:);
        v13 = node(Robin(:,3),:)-node(Robin(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        center = (node(Robin(:,1),:) + node(Robin(:,2),:)+node(Robin(:,3),:))/3;
        ii = find(isRobin);
        ss = pde.g_R(center).*area; % exact for linear g_R
        A = A + sparse(ii,ii,ss,Ndof,Ndof);
    end
    
    % Dirichlet boundary nodes: fixedFace
    fixedFace = []; freeFace= [];
    if ~isempty(bdFlag) % find boundary faces
        idxD = (bdFlag(:)==1); % all Dirichlet faces in bdFlag
        isFixedFace = false(NF,1);
        isFixedFace(elem2face(idxD)) = true;
        fixedFace = find(isFixedFace);
        freeFace = find(~isFixedFace);
    end
    
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given in the input
        s = accumarray(elem2face(:), 1, [NF 1]);
        fixedFace = find(s == 1);
        freeFace = find(s == 2);
    end 
    
    isPureNeumann = false;    
    if isempty(fixedFace) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedFace = 1;
        freeFace = 2:Ndof;    % eliminate the kernel by enforcing u(1) = 0;
    end
    
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedFace,fixedFace)=I, AD(fixedFace,freeFace)=0, AD(freeFace,fixedFace)=0.
    if ~isempty(fixedFace)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedFace) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
    %% Part 2: Find boundary faces and modify the right hand side b
    % Find boundary faces: Neumann
    Neumann = [];
    if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
        idxN = (bdFlag(:) == 2);      % all Neumann faces in bdFlag
        isNeumann = elem2face(idxN | idxR); % index of Neumann faces
        Neumann = face(isNeumann,:);      % Neumann faces        
    end
    
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag is given, only pde.g_N or pde.g_R is given in the input
        s = accumarray(elem2face(:), 1, [NF 1]);
        isNeumann = (s==1);         % index of boundary faces
        Neumann = face(s == 1,:);   % boundary faces are Neumann faces
    end
    
      % Neumann boundary condition
    if ~isempty(Neumann) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 3;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        gf = zeros(size(Neumann,1),1);
        v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
        v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:)...
                 + lambdagN(pp,3)*node(Neumann(:,3),:);
            gNp = pde.g_N(ppxy);
            gf = gf+ weightgN(pp)*gNp;
        end
        gf = gf.*area;
        b(isNeumann) = b(isNeumann) +gf;
    end
        
    % The case with non-empty Neumann faces but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann faces and no modification of
    % A,u,b is needed.
   
    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedFace) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        if isnumeric(pde.g_D)
            u(fixedFace) = pde.g_D(fixedFace);
        else % pde.g_D is a function handle
            center = (node(face(fixedFace,1),:)+node(face(fixedFace,2),:)+node(face(fixedFace,3),:))/3;
            u(fixedFace) = pde.g_D(center);
        end
        b = b - A*u;
        b(fixedFace) = u(fixedFace);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b); % compatilbe condition: sum(b) = 0
        b(1) = 0;
    end
    end % end of getbd3CR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of Poisson3CR