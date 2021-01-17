%% RATE OF CONVERGENCE OF PRIMAL-DUAL WEAK-GALERKIN for DIV-CURL SYSTEM
%
% This example is to show the rate of convergence on an L-shaped domain with a singular solution.
%
% See also DivCurl3PDWG
%
% Reference: A New Numerical Method for Div-Curl Systems with Low Regularity Assumptions
% S. Cao, C. Wang, J. Wang, https://arxiv.org/abs/2101.03466
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear; close all;
%%
maxIt = 4;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2UTotal = zeros(maxIt,1);
errqlambdaTotal = zeros(maxIt,1);
errSTotal = zeros(maxIt,1);
errL2QuTotal = zeros(maxIt,1); % \|Q_h u - u \|


%% Generate an initial mesh 
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1/2);
% [node,elem] = delmesh(node,elem,'(x < 0 & y > 0) | (z<0.5)');

% [node,elem] = delmesh(node,elem,'(x < 0 & y < 0) | (z<0) | (z>0.5)'); %
% for gradient

[node,elem] = delmesh(node,elem,'(x > 0 & y < 0) | (z<0) | (z>0.5)');
% for curl
bdFlag = setboundary3(node,elem,'Dirichlet');

%% data
pde = DataDivCurlSingular4;
option.inverseh = true;

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    [soln,eqn,info] = DivCurl3PDWG(node,elem,bdFlag,pde,option);
    N(k) = size(elem,1);
%     h(k) = (36./N(k)).^(1/3);
    h(k) = 1/2^k;
    
    errL2UTotal(k) = eqn.errorU;
    errqlambdaTotal(k) = eqn.errorqlambda;
    errSTotal(k) = eqn.errorS;
    
    [errL2QuTotal(k), ~] = getL2error3vec(node,elem,soln.u,pde);
%     [errL2QuTotal(k), ~] = getL2error3(node,elem,pde.exactu, soln.u);
   
    
    fprintf('PDWG method with h:     %3g \n', h(k));
    fprintf('Number of Dof:          %d \n', sum(eqn.freeDof));
    fprintf('L2 error in u:          %6g \n', errL2UTotal(k));
    fprintf('L2 error in Qu-u_h:     %6g \n', errL2QuTotal(k));
    fprintf('L2 error in (q,lambda): %6g \n', errqlambdaTotal(k));
    fprintf('L2 error in s:          %6g \n\n\n', errSTotal(k));
    
    % refine mesh
   if k < maxIt; [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag); end 
end

%%
visualizationPDwgDivCurl;