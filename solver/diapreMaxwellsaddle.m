function [u,p,info] =  diapreMaxwellsaddle(A,G,f,g,node,elem,HB,edge,isBdEdge,option)
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
D = diag(A);
B = tril(A);
Bt = B';
% Laplace for P1 element
Ap = assemblematrix(node,elem); % Lap for P1 element
setupOption.solver = 'NO';
setupOption.setupflag = 1;
isFreeNode = true(N,1);
isFreeEdge = ~isBdEdge;
isFreeNode(edge(isBdEdge,:)) = false;
setupOption.freeDof = isFreeNode;
Ap = Ap(isFreeNode,isFreeNode);
% setupOption.isFreeEdge = option.isFreeEdge;
if isempty(HB) 
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Ap,g,elem,setupOption);
else % HB is only needed for adaptive grid in 3D
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Ap,g,elem,setupOption,HB);
end
% transfer matrix
if Nu <= NE        % lowest order edge element
    II = node2edgematrix(node,edge,isBdEdge);
    II = II(isFreeEdge,repmat(isFreeNode,dim,1));
elseif Nu <= 2*NE  % first or second order edge element
    II = node2edgematrix1(node,edge,isBdEdge);
end
IIt = II';
grad = gradmatrix(edge,isBdEdge);
grad = grad(isFreeEdge,isFreeNode);
gradt = grad';

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
   option.solvermaxIt = 1; 
end
option.setupflag = false;
option.printlevel = 0;
option.x0 = zeros(Np,1);
% minres for the saddle point system
[x,flag,stopErr,itStep,err] = minres(bigA,bigF,tol,maxIt,@diagpreconditioner,[],x0);
% extract solution
u = x(1:Nu); 
p = x(Nu+1:end);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Diagonal Preconditioned MINRES \n');
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
        rc = reshape(IIt*ru(1:size(IIt,2)),Np,dim);  
        eaux = mg(Ap,rc,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        % transfer back to the edge element space
        eaux = [II*reshape(eaux,dim*Np,1); zeros(Nu-size(IIt,2),1)]; 
        % another correction from P1
        auxrp = gradt*ru;
        eaux2 = zeros(Nu,1);
        auxep = mg(Ap,auxrp,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        eaux2(1:size(grad,1)) = grad*auxep;
        e1 = eh + eaux + eaux2;
        % Laplacian for p
        e2 = mg(Ap,rp,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        e  =  [e1; e2];   
    end
end