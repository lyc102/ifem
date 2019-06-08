%% LSHAPE Problem
%
% LSHAPE solves the Poisson equation in a Lshaped domain.
%
% Detailed explanation can be found in 
%  <a href="matlab:ifem afemdoc">ifem afemdoc</a>
%
% See also  crack, Kellogg
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
%% Parameters
maxN = 5e4;     theta = 0.4;    maxIt = 100; 
N = zeros(maxIt,1);   
errL2 = zeros(maxIt,1);   
errH1 = zeros(maxIt,1); 
erreta = zeros(maxIt,1);

%%  Generate an initial mesh and set up PDE data
[node,elem] = squaremesh([-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');
pde = Lshapedata;

%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    [soln,eqn,info] = Poisson(node,elem,bdFlag,pde);
    figure(1);  showresult(node,elem,soln.u,[-50,12]);    
    % Step 2: ESTIMATE
%     eta = estimaterecovery(node,elem,u);         % recovery type
    eta = estimateresidual(node,elem,soln.u,pde);    % residual type
    % Record error and number of vertices
    errH1(k) = getH1error(node,elem,pde.Du,soln.Du);
    errL2(k) = getL2error(node,elem,pde.exactu,soln.u);
    erreta(k) = sqrt(sum(eta.^2));
    N(k) = size(node,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
end

%% Plot convergence rates
figure;
showrate3(N(1:k),errH1(1:k),10,'-*','||Du-Du_h||',...
          N(1:k),errL2(1:k),10,'k-+','||u-u_h||',...
          N(1:k),erreta(1:k),10,'m-x','\eta');