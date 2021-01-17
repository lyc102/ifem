%% RATE OF CONVERGENCE OF PRIMAL-DUAL WEAK-GALERKIN for DIV-CURL SYSTEM
%
% This example is to show the rate of convergence on a toroidal domain with 1 hole (b_1 = 1)
%
% See also DivCurl3PDWG
%
% Reference: A New Numerical Method for Div-Curl Systems with Low Regularity Assumptions
% S. Cao, C. Wang, J. Wang, https://arxiv.org/abs/2101.03466
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
clear; 
%%
maxIt = 3;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2UTotal = zeros(maxIt,1);
errqlambdaTotal = zeros(maxIt,1);
errL2QuTotal = zeros(maxIt,1); % \|Q_h u - u \|
errSTotal = zeros(maxIt,1);


%% Generate an initial mesh for 1 hole
[node,elem] = cubemesh([-2,2,-2,2,-2,2],1);
[node,elem] = delmesh(node,elem,...
    '(abs(x+0.5) <= 0.5 & abs(y+0.5)<=0.5) | (z<0) |(z>1) | (x>1) | (y>1)');
bdFlag = setboundary3(node,elem,'Dirichlet');
node = node/2;

%% data
gamma = 2/3; % 5/4, 1, 2/3
alpha = 2;
pde = DataDivCurlSingularTorus(gamma,alpha,2);
% pde = DataDivCurlLinear(1,1);
% pde = DataDivCurlSmooth(1,1);

%% Finite Element Method        
for k = 1:maxIt
    [soln,eqn,info] = DivCurl3PDWG(node,elem,bdFlag,pde);
    N(k) = size(elem,1);
    h(k) = (6./N(k)).^(1/3);
    errL2UTotal(k) = eqn.errorU;
    errqlambdaTotal(k) = eqn.errorqlambda;
    errSTotal(k) = eqn.errorS;
    
    
    [errL2QuTotal(k), ~, Qu] = getL2error3vec(node,elem,soln.u,pde);
    fprintf('\nPDWG method with h:     %3g \n', h(k));
    fprintf('Number of Dof:          %d \n', sum(eqn.freeDof));
    fprintf('L2 error in u:          %6g \n', errL2UTotal(k));
    fprintf('L2 error in Qu-u_h:     %6g \n', errL2QuTotal(k));
    fprintf('L2 error in (q,lambda): %6g \n', errqlambdaTotal(k));
    fprintf('L2 error in s:          %6g \n\n\n', errSTotal(k));
    
    if k < maxIt; [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag); end
end

%%
harmonic_u = Qu - soln.uh2elem;
visualizationPDwgDivCurl;
