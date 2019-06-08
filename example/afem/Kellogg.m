%% KEOLLOGG Problem
%
% KELLOGG solves a diffusion equation with jump coefficients with AFEM.
%
% KELLOGG solves the problem within maxN number of vertices. The
% input argument theta is a parameter used in the marking step. 
%
% See also  crack, Lshape
%
% Created by Chen-Song Zhang. Modified by Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all; clear variables;
%% Parameters
maxN = 3e3;     theta = 0.1;    maxIt = 1000; 
N = zeros(maxIt,1);     
uIuhErr = zeros(maxIt,1);
errH1 = zeros(maxIt,1);

%%  Generate an initial mesh
[node,elem] = squaremesh([-1 1 -1 1], 0.25);
bdFlag = setboundary(node,elem,'Dirichlet');

%% Set up PDE data
pde = Kelloggdata;

%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE/COARSEN*
for k = 1:maxIt
    % Step 1: SOLVE
    [soln,eqn,info] = Poisson(node,elem,bdFlag,pde);
    % Plot mesh and solution
    figure(1);  showresult(node,elem,soln.u,[27,26]);    
    % Step 2: ESTIMATE
%    eta = estimaterecovery(node,elem,u);            % recovery type
    eta = estimateresidual(node,elem,soln.u,pde);    % residual type
    % Record error and number of vertices
    uI = pde.exactu(node);
    uIuhErr(k) = sqrt((uI-soln.u)'*eqn.A*(uI-soln.u));
    errH1(k) = getH1error(node,elem,@pde.Du,soln.Du,5);
    N(k) = size(node,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem,bdFlag,HB,tree] = bisect(node,elem,markedElem,bdFlag);
    % COARSEN
    eta = eleminterpolate(eta,tree);
    markedElem = mark(elem,eta,0.5*theta,'COARSEN');
    [node,elem,bdFlag] = coarsen(node,elem,markedElem,bdFlag);
end

%%  Plot convergent rates in energy norm
figure(2)
showrate2(N(1:k),uIuhErr(1:k),20,'-*','|| Du_I-Du_h||',...
          N(1:k),errH1(1:k),20,'k-.','|| Du-Du_h||');