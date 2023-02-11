function [x,info] = amgMaxwellinterface(A,b,node,edge,option)
%% AMGMAXWELL algebraic multigrid solver for Maxwell equations.
% 
% x = amgMaxwell(A,b,node,edge) attempts to solve the system of
% linear equations Ax = b using multigrid type solver. The linear system is
% obtained by either the first or second family of linear edge element
% discretization of the Maxwell equation; See <a href="matlab:ifem
% coarsendoc">doc Maxwell</a>.
%
% amgMaxwell is more algebraic than mgMaxwell but still requires geometric
% information node and edge. Grapha Laplacian of vertices are used as
% auxiliary Poisson operator and amg is used as Poisson solver.
%
% Input 
%   -  A: curl(alpha curl) + beta I
%   -  b: right hand side
%   -  node,edge: mesh information
%   -  options: extra structures
%
% By default, the HX preconditioned PCG is used which works well for
% symmetric positive definite matrices (e.g. arising from eddy current
% simulation). For symmetric indefinite matrices (e.g. arising from time
% harmonic Maxwell equation), set option.solver = 'minres' or 'bicg' or
% 'gmres' to try other Krylov method with HX preconditioner.
%
% See also mg, mgMaxwell, mgMaxwellsaddle
%
% Reference: 
% R. Hiptmair and J. Xu, Nodal Auxiliary Space Preconditioning in
% H(curl) and H(div) Spaces. SIAM J. Numer. Anal., 45(6):2483-2509, 2007.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

t = cputime;
%% Initial check
Ndof = length(b);                  % number of dof
N = size(node,1);                  % number of nodes
NE = size(edge,1);                 % number of edge;
dim = size(node,2);      
% Assign default values to unspecified parameters
if ~exist('option','var')
    option = []; 
end 
option = mgoptions(option,length(b));    % parameters
x0 = option.x0; 
% N0 = option.N0; 
tol = option.tol; 
%tol = 5*10^(-7);
% mu = option.smoothingstep; preconditioner = option.preconditioner; coarsegridsolver = option.coarsegridsolver; 
printlevel = option.printlevel; %setupflag = option.setupflag;
maxIt = 2000; %400;   % increase the default step (200) for Maxwell equations
% tol = 1e-7;    % reset the tol
% Check for zero right hand side
if (norm(b) == 0)                  % if rhs vector is all zeros
    x = zeros(Ndof,1);             % then solution is all zeros
    flag = 0;                      
    itStep = 0;                    
    err = 0;  
    time = cputime - t;
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));        
    return
end

%% Transfer operators from nodal element to edge element
if isfield(option,'isBdEdge')
    isBdEdge = option.isBdEdge;
else
    deg = sum(spones(A(1:NE,1:NE)),2);
    isBdEdge = (deg == 1);
end
if Ndof == NE        % lowest order edge element
    II = node2edgematrix(node,edge,isBdEdge);
elseif Ndof >= 2*NE  % first or second order edge element
    II = node2edgematrix1(node,edge,isBdEdge);
end
IIt = II';
[grad,isBdNode] = gradmatrix(edge,isBdEdge);
gradt = grad';

%% Block smoother
%if option.blklevel == 0 || ~isfield(option,'blklevel')
    
interfaceEdgeID = option.blkId+1:NE;
interfaceEdge = false(NE, 1);
interfaceEdge(option.blkId+1:NE) = true;
regularEdge = false(NE, 1);
regularEdge(1:option.blkId) = true;

if isfield(option,'blklevel')
    
    ExtLevel = 0;
    while (ExtLevel<option.blklevel)
        
        V2Emat= sparse([1:NE,1:NE]',[edge(:,1);edge(:,2)],1,NE,N);
        interfaceNodeID = union(edge(interfaceEdgeID,1),edge(interfaceEdgeID,2));
        interfaceNode = false(N, 1); interfaceNode(interfaceNodeID) = true;
        interfaceEdge = (V2Emat*interfaceNode>0);
        regularEdge = ~interfaceEdge;
        interfaceEdgeID = find(interfaceEdge==1);
        ExtLevel = ExtLevel+1;
        
    end
    
end


Ainterface = A(interfaceEdge,interfaceEdge);
Aregular = A(regularEdge,regularEdge);
B = A(interfaceEdge, regularEdge);
D = diag(Aregular);

tic
if isfield(option,'fact') == 1
    switch option.fact
        case 'lu'
            [L,U,P,Q] = lu(Ainterface);
        case 'chol'
            [R,flag,P] = chol(Ainterface);
    end
else
    [L,U,P,Q] = lu(Ainterface);
end
toc

%% Auxiliary Poisson matrix
%   -  A: curl(alpha curl) + beta I
%   - AP: - div(alpha grad) + beta I
%   - BP: - div(beta grad)

if isfield(option,'AP')
    AP = option.AP;
   edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
%     edgeLength = sqrt(sum(edgeVec.^2,2));
%     if isfield(option,'beta') % resacle by the dielectric coefficients
%         if isreal(option.beta) && (length(option.beta) == NE)
%            beta = option.beta;  
%         else % option.beta is a function 
%            edgeMiddle = (node(edge(:,2),:) + node(edge(:,1),:))/2; 
%            beta = option.beta(edgeMiddle);         
%         end
%     end
%     M = accumarray(edge(:),repmat((edgeLength.^3).*beta,2,1),[N 1]);
%     AP = AP + spdiags(M,0,N,N);
else
    % build graph Laplacian to approximate AP
    edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
    edgeLength = sqrt(sum(edgeVec.^2,2));
    if isfield(option,'alpha') % resacle by the magnetic coefficients
        if isreal(option.alpha) && (length(option.alpha) == NE)
           alpha = option.alpha;  
        else % option.alpha is a function 
           edgeMiddle = (node(edge(:,2),:) + node(edge(:,1),:))/2; 
           alpha = option.alpha(edgeMiddle);         
        end
    end
    % edge weight: h*alpha
    % AP = gradt*spdiags(edgeLength.*alpha,0,NE,NE)*grad;
    i1 = (1:NE)'; j1 = double(edge(:,1)); s1 = sqrt(edgeLength.*alpha);
    i2 = (1:NE)'; j2 = double(edge(:,2)); s2 = -s1;
    isFreeEdge = ~isBdEdge;
    G = sparse([i1(isFreeEdge);i2(isFreeEdge)],...
        [j1(isFreeEdge);j2(isFreeEdge)],...
        [s1(isFreeEdge);s2(isFreeEdge)],NE,N);
    %i1 = (1:NE)'; j1 = double(edge(:,1)); s1 = ones(size(alpha)).*sqrt(alpha);
%     i2 = (1:NE)'; j2 = double(edge(:,2)); s2 = -s1;
%     isFreeEdge = ~isBdEdge;
%     G = sparse([i1(:);i2(:)],...
%         [j1(:);j2(:)],...
%         [s1(:);s2(:)],NE,N);
    AP = G'*G;
    %E = G'*A*G;
    % lumped mass matrix: h^3*beta
    if isfield(option,'beta') % resacle by the dielectric coefficients
        if isreal(option.beta) && (length(option.beta) == NE)
           beta = option.beta;  
        else % option.beta is a function 
           edgeMiddle = (node(edge(:,2),:) + node(edge(:,1),:))/2; 
           beta = option.beta(edgeMiddle);         
        end
    end
    M = accumarray(edge(:),repmat((edgeLength.^3).*beta,2,1),[N 1]);
    AP = AP + spdiags(M,0,N,N);
end
% BP is Galerkin projection to the free node space
% boundary nodes
bdidx = zeros(N,1); 
bdidx(isBdNode) = 1;
Tbd = spdiags(bdidx,0,N,N);
BP = gradt*A(1:NE,1:NE)*grad + Tbd;
%BP = gradt*spdiags(edgeLength.*beta,0,NE,NE)*grad + Tbd;
% if strcmp(option.outsolver,'minres')
%     BP = gradt*option.BPP(1:NE,1:NE)*grad + Tbd;
% end

setupOption.solver = 'NO';
[x,info,APi,Ri,RRi,ResAP,ProAP,clA] = amg(AP,ones(N,1),setupOption); %#ok<ASGLU>
[x,info,BPi,Si,SSi,ResBP,ProBP,clB] = amg(BP,ones(N,1),setupOption); %#ok<ASGLU>

%% Krylov iterative methods with HX preconditioner
k = 1;
err = 1;
switch upper(option.outsolver)
    case 'CG'
        if printlevel>=1
            fprintf('Conjugate Gradient Method using HX preconditioner \n');
        end
        x = x0;
        r = b - A*x;
        nb = norm(b);
        err = zeros(maxIt,2);
        err(1,:) = norm(r)/nb; 
        while (max(err(k,:)) > tol) && (k <= maxIt)
            % compute Br by HX preconditioner
            Br = HXpreconditioner(r); 
            % update tau, beta, and p
            rho = r'*Br;  % r'*Br = e'*ABA*e approximates e'*A*e
            if k==1
                p = Br;
            else
                beta = rho/rho_old;
                p = Br + beta*p;
            end
            % update alpha, x, and r
            Ap = A*p;
            alpha = rho/(Ap'*p);
            r = r - alpha*Ap;
            x = x + alpha*p;
            rho_old = rho;
            % compute err for the stopping criterion
            k = k + 1;
            err(k,1) = sqrt(abs(rho/(x'*b))); % approximate relative error in energy norm
            err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
            if printlevel >= 2
                fprintf('#dof: %8.0u, HXCG iter: %2.0u, err = %12.8g\n',...
                         Ndof, k, max(err(k,:)));
            end
        end
        err = err(1:k,:);
        itStep = k-1;
        if k > maxIt || (max(err(end,:))>tol)
            flag = 1;
        else
            flag = 0;
        end
    case 'MINRES'
        fprintf('Minimum Residual Method with HX preconditioner \n')
        [x,flag,err,itStep] = minres(A,b,tol,maxIt,@HXpreconditioner,[],x0);         
%         x = minres(A,b,tol,maxIt,@HXpreconditioner,[],x0);         
    case 'GMRES'
        fprintf('General Minimum Residual Method with HX preconditioner \n')
        tic
        [x,flag,err,itStep] = gmres(A,b,10,tol,maxIt,@HXpreconditioner,[],x0);
        toc
        itStep = (itStep(1) -1)*10 + itStep(2);
end

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('#dof: %8.0u,   #nnz: %8.0u,   iter: %2.0u,   err = %8.4e,   time = %4.2g s\n',...
                 Ndof, nnz(A), itStep, max(err(end,:)), time)
end
if  (flag == 1) && (printlevel>0)
    fprintf('NOTE: the iterative method does not converge');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions HXpreconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Br = HXpreconditioner(r)
    %% 1. Smoothing in the finest grid of the original system
    if strcmp(option.smoother,'BD')
        if isfield(option,'fact') == 1
            switch option.fact
                case 'lu'
                    eh = zeros(NE, 1);
                    eh(regularEdge,:) = 0.75*r(regularEdge,:)./D;
                    %        eh(interfaceEdge,:) = Ainterface\r(interfaceEdge,:);
                    rtmp = P*r(interfaceEdge,:);
                    rtmp = L\rtmp;
                    rtmp = U\rtmp;
                    rtmp = Q*rtmp;
                    eh(interfaceEdge,:) = rtmp;
                case 'chol'
                    eh = zeros(NE, 1);
                    eh(regularEdge,:) = 0.75*r(regularEdge,:)./D;
                    %        eh(interfaceEdge,:) = Ainterface\r(interfaceEdge,:);
                    rtmp = P'*r(interfaceEdge,:);
                    rtmp = R'\rtmp;
                    rtmp = R\rtmp;
                    rtmp = P*rtmp;
                    eh(interfaceEdge,:) = rtmp;
            end
        else
            eh = zeros(NE, 1);
            eh(regularEdge,:) = 0.75*r(regularEdge,:)./D;
            %        eh(interfaceEdge,:) = Ainterface\r(interfaceEdge,:);
            rtmp = P*r(interfaceEdge,:);
            rtmp = L\rtmp;
            rtmp = U\rtmp;
            rtmp = Q*rtmp;
            eh(interfaceEdge,:) = rtmp;
        end
    elseif strcmp(option.smoother,'GS')
        %         eh = triu(A)\(D.*(tril(A)\r)); # does not work if interface edge is
        %         not treated differently
        rN = r(regularEdge); rI = r(interfaceEdge);
        rN = rN./D; rI = Ainterface\rI;
        ehN = rN + B'*(Ainterface\(B*rN))./D - B'*rI./D;
        ehI = -Ainterface\(B*rN) + rI;
        eh = [ehN; ehI];
    else
        eh = 0.75*r./D;  % Jacobi method. less computational time
    end
    %% 2. amg solver for auxiliary operators
    amgoption.solvermaxit = 3;  % 3 for SPD matrix
    amgoption.printlevel = 0;
    rc = reshape(IIt*r(1:size(IIt,2)),N,dim);   % transfer to the nodal linear element space
    %eaux = II*reshape(AP\rc,dim*N,1);
    %eaux = II*reshape(amg(AP,rc,amgoption),dim*N,1);
    rAP = amgInterface(AP,rc,APi,Ri,RRi,ResAP,ProAP,clA,amgoption);
    eaux = II*reshape(rAP,dim*N,1);
    eh = eh + eaux;
    rb = gradt*r(1:NE);
    %eauxb = grad*(BP\rb);
    %eauxb = grad*(amg(BP,rb,amgoption));
    rBP = amgInterface(BP,rb,BPi,Si,SSi,ResBP,ProBP,clB,amgoption);
    eauxb = grad*rBP;
    Br = eh + eauxb;
    end

    function z = amgInterface(M,r0,Ai,Bi,BBi,Res,Pro,cl,amgoption)
        
        maxItz = amgoption.solvermaxit;
        Nb = size(r0,2);
        z = repmat(zeros(size(r0,1),1),1,Nb);
        kz = 1;
        rz = r0 - M*z;
        nbz = max(sqrt(sum(r0.^2,1)));
        errz = zeros(maxItz,2);
        errz(1,:) = max(sqrt(sum(rz.^2,1)))/nbz;

        amgprintlevel = amgoption.printlevel;  
        level = max(min(ceil(log2(N)/2-4),8),2);
        prefunc = @wcycle;
        if amgprintlevel >= 1
            fprintf('Conjugate Gradient Method\n')
        end
        if isfield(amgoption,'smoothingstep')  % smoothing steps
            mu = amgoption.smoothingstep;
        else
            mu = 2;
        end
      
        
        while (max(errz(kz,:)) > tol) && (kz <= maxItz)    
            % compute Br by MG
            Brz = prefunc(r0);
            % update tau, beta, and p
            rhoz = dot(Brz,r0);  % e'*ABA*e approximates e'*A*e
            if kz == 1
                pz = Brz;
            else
                betaz = rhoz./rho_oldz;
                pz = Brz + betaz.*pz;
            end
            % update alpha, x, and r
            Apz = M*pz;
            alpha = rhoz./dot(Apz,pz);
            r0 = r0 - alpha.*Apz;
            z = z + alpha.*pz;
            rho_oldz = rhoz;
            kz = kz + 1;
            % compute err for the stopping criterion
        %     err(k,1) = alpha*sqrt(p'*Ap/(x'*A*x)); % increamental error in energy norm
            errz(kz,1) = max(sqrt(abs(rhoz./dot(z,r0)))); % approximate relative error in energy norm
            errz(kz,2) = max(sqrt(sum(r0.^2,1)))/nbz; % relative error of the residual in L2-norm
            if amgprintlevel >= 2
                fprintf('#dof: %8.0u, MGCG iter: %2.0u, err = %8.4e\n',...
                         N, kz-1, max(errz(kz,:)));
            end
        end
        errz = errz(1:kz,:);
        itStepz = kz-1;
        
        function e = wcycle(r,J)        % solve equations Ae = r in each level
            if nargin<=1
                J = level;
            end
            if J == cl
                e = Ai{cl}\r;   % exact solver in the coaresest grid
                return
            end
            % fine grid pre-smoothing
            e = Bi{J}\r;   % pre-smoothing
            for s = 1:mu           % extra mu steps smoothing
                e = e + Bi{J}\(r-Ai{J}*e);
            end
            rc = Res{J}*(r - Ai{J}*e);
            % coarse grid correction twice
            ec = wcycle(rc,J-1);
            ec = ec + wcycle(rc - Ai{J-1}*ec,J-1);
            % fine grid post-smoothing
            e = e + Pro{J-1}*ec;
            e = e + BBi{J}\(r-Ai{J}*e);
            for s = 1:mu
                e = e + BBi{J}\(r-Ai{J}*e); % post-smoothing
            end
        end
        
    end

end