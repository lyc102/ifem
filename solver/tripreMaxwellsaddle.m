function [u,p,info] =  tripreMaxwellsaddle(A,G,f,g,node,elem,bdFlag,M,grad,option)
%Solve the maxwell system with divgence free condition,
%         [A  G] [u]  = f                (1.1)
%         [G' O] [p]  = g0               (1.2)
% where  G = M_e*grad.
% This system can be rewritten as 
%         [A+G*DMinv*G'       G] [u]  = f +G*DMinv*g0
%         [G'                 O] [p]  = g0
% in fact that
%  [A+G*DMinv*G'  G] [I    grad]                = [Abar   O  ]
%  [G'            O] [0  -Dminv*grad'*M_e*grad ]= [G'     A_p]
%
%  We solve the (1.1)-(1.2) with the GMRES, the preconditioner is
%
%       [I   grad                ]  [Abar O   ]^{-1}
%       [O  -Dminv*grad'*M_e*grad]  [G'   A_p ]
%
% Created by Long chen and Jie Zhou on Aug,2015.


Nf = length(f); 
Ng = length(g);
t = cputime;

%% Parameters
if ~exist('option','var'), option = []; end
Ndof = Nf + Ng; 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; 
tol = option.tol; 
maxIt = option.solvermaxit; 
printlevel = option.printlevel; 
d = size(node,2);

%% Set up auxiliary matrices
if d == 2
    area = simplexvolume(node,elem);
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,...
                        [max(elem(:)),1]);
elseif d == 3
    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation presereved
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);    
end
DMinv = spdiags(1./Mvlump(option.isFreeNode),0,Ng,Ng);
f = f + G*(DMinv*g);  % add second equation to the first one
Abar = A + G*DMinv*G'; % Hodge Laplacian
Ap = grad'*M*grad;    % scalar Laplacian

%% Set up matrices for multigrid
setupOption.solver = 'NO';
setupOption.isFreeEdge = option.isFreeEdge;
setupOption.freeDof = option.isFreeNode;
[~,~,Ai_N,Bi_N,BBi_N,Res_N,Pro_N] = mg(Ap,g,elem,setupOption);
[x,info,Ai,Bi,BBi,Res,Pro] = mgHodgeLapE(Abar,f,node,elem,bdFlag,setupOption); %#ok<*ASGLU>

%% Form a big matrix equation
bigA = [Abar G; G' sparse(Ng,Ng)];
bigF = [f; g];

%% Preconditioned GMRES
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
option.x0 = zeros(Nf,1);
option.printlevel = 0;
Apmgoption.printlevel = 0;
Apmgoption.printlevel = 0;
Apmgoption.setupflag = false;
Apmgoption.maxIt     = 1;
% minres for the saddle point system
[x,flag,stopErr,itStep,err] = gmres(bigA,bigF,20,tol,maxIt,@tripre,[],x0);
itStep = (itStep(1)-1)*20 + itStep(2); % total iteration
u = x(1:Nf);
p = x(Nf+1:end);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Triangular Preconditioned GMRES \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(bigA), option.solvermaxIt, itStep, stopErr, time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',stopErr);
  
%% Preconditioner   
    function e  = tripre(r)
        r1 = r(1:Nf);
        r2 = r(Nf+1:end);
        e1 = mgHodgeLapE(Abar,r1,node,elem,bdFlag,option,Ai,Bi,BBi,Res,Pro); 
        e2 = mg(Ap,r2-G'*e1,elem,Apmgoption,Ai_N,Bi_N,BBi_N,Res_N,Pro_N);
        e  = [e1 + grad*e2; -DMinv*(Ap*e2)];
    end
end