function [u,p,info] =  diapreMaxwellsaddle(A,G,f,g,node,elem,edge,option,varargin)
%Solve the maxwell system with divgence free condition,
%         [A  G] [u]  = f                (1.1)
%         [G' O] [p]  = g                (1.2)
% where  G = M_e*grad.
%



Nf = length(f); 
Ng = length(g);
N = size(node,1);
dim = size(node,2);
t = cputime;

%% Parameters
if ~exist('option','var'), option = []; end
Ndof = Nf + Ng; 
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
% setupOption.isFreeEdge = option.isFreeEdge;
HB = varargin{end};
if isempty(HB) 
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Ap,g,elem,setupOption);
else % HB is only needed for adaptive grid in 3D
    [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(Ap,g,elem,setupOption,HB);
end
% transfer matrix
if Ndof == NE        % lowest order edge element
    II = node2edgematrix(node,edge,isBdEdge);
elseif Ndof >= 2*NE  % first or second order edge element
    II = node2edgematrix1(node,edge,isBdEdge);
end
IIt = II';

%% Form a big matrix equation
bigA = [A,G; G',sparse(Ng,Ng)];
bigF = [f; g];

%% Preconditioned MINRES
% options for V-cycle of Schur complement
if strcmp(option.solver,'CG') % change default set up in mgoption 
   option.solver = 'Vcycle';
end
if isfield(option,'Vit')
   option.solvermaxIt = option.Vit;
else
   option.solvermaxIt = 1; 
end
option.setupflag = false;
option.x0 = zeros(Ndof,1);
option.printlevel = 0;
% minres for the saddle point system
[x,flag,stopErr,itStep,err] = minres(bigA,bigF,tol,maxIt,@diagpreconditioner,[],x0);
% extract solution
u = x(1:Nf); 
p = x(Nf+1:end);

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
        ru = r(1:Nf);
        rp = r(Nf+1:end);  
        % HX preconditioner for curlcurl matrix
        eh = Bt\(D.*(B\ru));  % Gauss-Seidal for A
        % transfer to the vector P1 linear element space
        rc = reshape(IIt*ru(1:size(IIt,2)),N,dim);  
        eaux = mg(Ap,rc,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        % transfer back to the edge element space
        eaux = [II*reshape(eaux,dim*N,1); zeros(Nf-size(IIt,2),1)]; 
        e1 = eh + eaux;
        % Laplacian for p
        e2 = mg(Ap,rp,elem,option,Ai,Bi,BBi,Res,Pro,isFreeDof);
        e  =  [e1; e2];   
    end
end