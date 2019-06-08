function [sigma,u,info] =  sdiagpreHodgeLapE(M,G,C,f,g,node,elem,bdFlag,option)
%% SDIAPREHODGELAPE: solve the Hodge Laplacian by a signed diagonal preconditioner 
%
%
%
% Rewritten by Long Chen on Feb 21, 2014.


t = cputime;
%% Parameters
if ~exist('option','var'), option = []; end
Ndof = length(f)+length(g); 
option = mgoptions(option,Ndof);    % parameters
x0 = option.x0; 
tol = option.tol; 
maxIt = option.solvermaxit; 
printlevel = option.printlevel; 
d = size(node,2);

%% Set up auxiliary matrices
Nf = length(f); Ng = length(g);
if d == 2
    area = simplexvolume(node,elem);
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,...
                        [max(elem(:)),1]);
elseif d == 3
    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation presereved
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);    
end
DMinv = spdiags(1./Mvlump(option.isFreeNode),0,Nf,Nf);
Gt = G';
Abar = G*DMinv*Gt + C;

%% Set up matrices for multigrid
setupOption.solver = 'NO';
setupOption.isFreeEdge = option.isFreeEdge;
[x,info,Ai,Bi,BBi,Res,Pro] = mgHodgeLapE(Abar,g,node,elem,bdFlag,setupOption); %#ok<*ASGLU>

%% Form a big matrix equation
bigA = [-M,Gt; G,C];
bigF = [f; g+G*(DMinv*f)];
% modify bigA to L*bigA*U
% B = G*(speye(size(DMinv)) - DMinv*M);
% bigA = [M, -B'; B, B*DMinv*Gt + Abar];
% bigF = [-f; g+G*(DMinv*f)];

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
option.x0 = zeros(Ng,1); % set to zero for one V-cycle preconditioner
option.printlevel = 0;   % x0 is used in gmres input
% minres for the saddle point system
[x,flag,stopErr,itStep,err] = gmres(@LbigAU,bigF,20,tol,maxIt,@sdiagpre,[],x0);
itStep = (itStep(1)-1)*20 + itStep(2); % total iteration
sigma = x(1:Nf);
u = x(Nf+1:end);
sigma = sigma + DMinv*(Gt*u);

%% Output
time = cputime - t;
if printlevel >= 1
    fprintf('Signed Diagonal Preconditioned GMRES \n');
    fprintf('#dof: %8.0u,  #nnz: %8.0u, V-cycle: %2.0u, iter: %2.0u,   err = %8.2e,   time = %4.2g s\n',...
                 Ndof, nnz(bigA), option.solvermaxIt, itStep, stopErr, time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',stopErr);

%% Modifed matrix
    function b = LbigAU(x)
        x(1:Nf) = x(1:Nf) + DMinv*(Gt*x(Nf+1:end));
        b = bigA*x;
        b(Nf+1:end) = G*(DMinv*b(1:Nf)) + b(Nf+1:end);
    end

%% Preconditioner 
    function e = sdiagpre(r)
        r1 = r(1:Nf);
        r2 = r(Nf+1:end);
%         r2 = G*(DMinv*r1) + r2;
        e1 = -DMinv*r1;
        e2 = mgHodgeLapE(Abar,r2,node,elem,bdFlag,option,Ai,Bi,BBi,Res,Pro);
        e  = [e1; e2];
    end
end