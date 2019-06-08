function [u,p,info] =  diapreMaxwellsaddle(A,G,f,g,node,elem,bdFlag,option)
%Solve the maxwell system with divgence free condition,
%         [A  G] [u]  = f                (1.1)
%         [G' O] [p]  = g0               (1.2)
% where  G = M_e*grad.
%
% This system can be rewritten as 
%         [A+G*DMinv*G'       G] [u]  = f +G*DMinv*g0
%         [G'                 O] [p]  = g0
%
%  We solve the above system using minres with the diagonal preconditioner
%
%      [Abar O]^{-1}
%      [ O   M]
%
% Created by Jie Zhou based on tripreMaxwellsaddle on Aug 23, 2015.
%
% Modified by Long Chen. Sep 29, 2015. Use a different diagonal
% precondtioner motivated by Stokes equation.


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
    volume = abs(simplexvolume(node,elem)); % uniform refinement in 3D is not orientation presereved
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);    
end
DMinv = spdiags(1./Mvlump(option.isFreeNode),0,Ng,Ng);
% f = f + G*(DMinv*g);  % add second equation to the first one
Abar = G*DMinv*G' + A;

%% Set up matrices for multigrid
setupOption.solver = 'NO';
setupOption.isFreeEdge = option.isFreeEdge;
[x,info,Ai,Bi,BBi,Res,Pro] = mgHodgeLapE(Abar,f,node,elem,bdFlag,setupOption); %#ok<*ASGLU>

%% Form a big matrix equation
% bigA = [Abar,G; G',sparse(Ng,Ng)];
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
option.x0 = zeros(Nf,1);
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
        r1 = r(1:Nf);
        r2 = r(Nf+1:end);          
        e1 =  mgHodgeLapE(Abar,r1,node,elem,bdFlag,option,Ai,Bi,BBi,Res,Pro);
        e2 =  DMinv*r2;              
        e  =  [e1; e2];   
    end
end