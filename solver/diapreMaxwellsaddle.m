function [u,p,info] =  diapreMaxwellsaddle(A,G,f,g,node,elem,Me,HB,edge,isBdEdge,option)
%Solve the maxwell system with divgence free condition,
%         [A  G] [u]  = f                (1.1)
%         [G' O] [p]  = g                (1.2)
% where  G = M_e*grad.
%



t = cputime;

%% Parameters and options
Nu = length(f); 
Np = length(g);
N = size(node,1);
NE = size(edge,1);
dim = size(node,2);
Ndof = Nu + Np; 
if ~exist('option','var'), option = []; end
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; 
tol = option.tol; 
maxIt = option.solvermaxit; 
printlevel = option.printlevel; 
dim = size(node,2);

%% Set up auxiliary matrices
AM = A + Me;
D = diag(AM);
B = tril(AM);
Bt = B';
% Laplace for P1 element or P2 element
Ap1 = assemblematrix(node,elem); % Lap for P1 element
if  Np > N && Np <= N+NE % Lap for P2 element
    mesh.node = node;
    mesh.elem = elem;
    mesh.size = size(node,1);
    if dim == 2
        eqn = getfemmatrix(mesh,'Poisson','P2');
        Ap = eqn.A;
    end
    if dim == 3
        eqn = getfemmatrix3(mesh,'Poisson','P2');
        Ap = eqn.A;
    end
end
setupOption.solver = 'NO';
setupOption.setupflag = 1;
isFreeNode = true(N,1);
isFreeEdge = ~isBdEdge;
isFreeNode(edge(isBdEdge,:)) = false;
Ap1 = Ap1(isFreeNode,isFreeNode);
Np1 = size(Ap1,1);  % length of free nodes
Ne1 = sum(isFreeEdge); % length of free edges
setupOption.freeDof = isFreeNode;
if Np > N && Np <= N+NE  % P2 element
    freeP2Dof = [isFreeNode; isFreeEdge];
    Ap = Ap(freeP2Dof,freeP2Dof);
end
% setupOption.isFreeEdge = option.isFreeEdge;
if isempty(HB) 
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Ap1,zeros(Np1,1),elem,setupOption);
else % HB is only needed for adaptive grid in 3D
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Ap1,zeros(Np1,1),elem,setupOption,HB);
end
% transfer matrix
if Nu <= NE && Np <= N   % lowest order edge element
    II = node2edgematrix(node,edge,isBdEdge);
    II = II(isFreeEdge,repmat(isFreeNode,dim,1));
elseif Nu <= 2*NE && Np <= N+NE  % first or second order edge element
    II = node2edgematrix1(node,edge,isBdEdge);
    II = II([isFreeEdge; isFreeEdge],repmat(isFreeNode,dim,1));
end
IIt = II';
grad = gradmatrix(edge,isBdEdge);
grad = grad(isFreeEdge,isFreeNode);
gradt = grad';
if Np > N && Np <= N+NE  % P2 element
    P1toP2 = sparse([(1:N)'; N+(1:NE)'; N+(1:NE)'], ...
                    [(1:N)'; double(edge(:))],...
                    [ones(N,1); 0.5*ones(2*NE,1)]',Ndof,N);
    auxPro = P1toP2(freeP2Dof,isFreeNode);        
    auxRes = auxPro';
end

%% Form a big matrix equation
bigA = [A,G; G',sparse(Np,Np)];
bigF = [f; g];

%% Preconditioned MINRES
% options for V-cycle of Ap
if strcmp(option.solver,'CG') % change default set up in mgoption 
   option.solver = 'Vcycle';
end
if isfield(option,'Vit')
   option.solvermaxIt = option.Vit;
else
   option.solvermaxIt = 2; 
end
option.setupflag = false;
option.printlevel = 0;
option.x0 = zeros(Np1,1);
% minres for the saddle point system
[x,flag,stopErr,itStep,err] = minres(bigA,bigF,tol,maxIt,@diagpreconditioner,[],x0);
% extract solution
u = x(1:Nu); 
p = x(Nu+1:end);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Diagonal Preconditioned MINRES for Maxwell Equations \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(bigA), option.solvermaxIt, itStep, stopErr, time)
end

if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',stopErr);

%% Preconditioner 
    function e = diagpreconditioner(r)
        % separate residual 
        ru = r(1:Nu);
        rp = r(Nu+1:end);  
        % HX preconditioner for curlcurl matrix
        eh = Bt\(D.*(B\ru));  % Gauss-Seidal for A
        % transfer to the vector P1 linear element space
        rc = reshape(IIt*ru,Np1,dim);  
        eaux = mg(Ap1,rc,elem,option,Ai,Bi,BBi,Res,Pro);
        % transfer back to the edge element space
        eaux = II*reshape(eaux,dim*Np1,1); 
        % another correction from P1
        auxrp = gradt*ru(1:Ne1);
        eaux2 = zeros(Nu,1);
        auxep = mg(Ap1,auxrp,elem,option,Ai,Bi,BBi,Res,Pro);
        eaux2(1:Ne1) = grad*auxep;
        e1 = eh + eaux + eaux2;
        % Laplacian for p
        if Np <= N
            e2 = mg(Ap1,rp,elem,option,Ai,Bi,BBi,Res,Pro);
        elseif Np <= N+NE
            e2 = triu(Ap)\(diag(Ap).*(tril(Ap)\rp));
            rc = auxRes*(rp-Ap*e2);
            ep1 = mg(Ap1,rc,elem,option,Ai,Bi,BBi,Res,Pro);
            e2 = e2 + auxPro*ep1;
        end
        e  =  [e1; e2];   
    end
end