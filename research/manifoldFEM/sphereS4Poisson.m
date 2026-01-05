%% SPHERES4POISSON Poisson equation on S^4 using stereographic projection and tensorial elements
%   the quad-linear finite element approximation of the Poisson equation
%
%       -Delta u + b*u = f  on S^4
%
%   with no boundary conditions.
%
%   This solver is designed for solving subproblems on each chart (a hyper-cube domain) in the
%   domain decomposition of manifold PDEs using stereographic projection.
%   The boundary values are transferred from the other chart using
%   interpolation weights. 
%
% See also
%   cubePoissonQ1, PoissonQ1

% Reference:
%   Cao, Shuhao, and Lizhen Qin. "A numerical domain decomposition method for solving elliptic equations on manifolds."
%   SIAM Journal on Scientific Computing 46, no. 1 (2024): A376-A398.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

clear;
dim = 4; % S^4, can be embedded to R^5
b = 1; % L2 coeff
n = 16;
r = 1 + 0.5/n; % overlapping radius needs > 1 (this choice overlapping grid is one element)
tol = 1e-8;
maxIter = 1e4;

close all;
clc;

%% Setup mesh and node information
pde = S4LaplacianBeltramidata;
[node, isBdNode, isUpperBoundary] = hypercubemesh(dim, n);

elemNonOverlap = find(~isUpperBoundary); % idx of first node of each element that does not touch upper boundary
nodeInt = node(~isUpperBoundary,:);

h=2*r/n; % mesh size
center = -r + (nodeInt + 0.5) * h; % element center
vol = getVolS4(dim,center);
A = assembleStiffness(pde,center,elemNonOverlap,h,n,vol);

rhsUpper = assembleRHS(pde, center, elemNonOverlap, h, n, vol, 'upper');
rhsLower = assembleRHS(pde, center, elemNonOverlap, h, n, vol, 'lower');
%%
bdNode=find(isBdNode);
intNode=find(~isBdNode);
[bdNodeIdxOther, bdInterpWeights] = transitionChartMap(node, bdNode, r, n);

uhUpper= zeros((n+1)^dim, 1);
uhLower= zeros((n+1)^dim, 1);

uIUpper = pde.exactu(node, r, h, 'upper');
uILower = pde.exactu(node, r, h, 'lower');

errUpper = uIUpper - uhUpper;
errLower = uILower - uhLower;
ResiduaLInftyErr = max(abs([errUpper',errLower']));
ResiduaLInftyErrDiff = 1e1;
%%
iter=0;
while ResiduaLInftyErrDiff > tol
    uhUpper(bdNode) =0;
    for j=1:2^dim
        % uhUpper(bdNode) = uhUpper(bdNode) + uhLower(BoundaryInf(:,2*j -1)).* BoundaryInf(:,2*j);
        uhUpper(bdNode) = uhUpper(bdNode) + uhLower(bdNodeIdxOther(:,j)).* bdInterpWeights(:,j);
    end
    rhsUpperRes = rhsUpper(intNode) - A(intNode, bdNode) * uhUpper(bdNode);

    if size(node,1) < 1e4
        uhUpper(intNode) = A(intNode, intNode) \ rhsUpperRes;
    else
        uhUpper(intNode) = pcg(A(intNode, intNode), rhsUpperRes,tol,maxIter,[],[],uhUpper(intNode));
    end


    uhLower(bdNode) =0;
    for j=1:2^dim
        % uhLower(bdNode) = uhLower(bdNode) + uhUpper(BoundaryInf(:,2*j -1)).* BoundaryInf(:,2*j);
        uhLower(bdNode) = uhLower(bdNode) + uhUpper(bdNodeIdxOther(:,j)).* bdInterpWeights(:,j);
    end
    rhsLowerRes = rhsLower(intNode) - A(intNode, bdNode) * uhLower(bdNode);

    if size(node,1) < 1e4
        uhLower(intNode) = A(intNode, intNode) \ rhsLowerRes;
    else
        uhLower(intNode) = pcg(A(intNode, intNode), rhsLowerRes, tol,maxIter,[],[],uhLower(intNode));
    end

    iter=iter+1;
    errUpper = uIUpper - uhUpper;
    errLower = uILower - uhLower;
    ResidualInftyErrPrev = ResiduaLInftyErr;
    ResiduaLInftyErr = max(abs([errUpper',errLower']));
    ResiduaLInftyErrDiff = abs(ResiduaLInftyErr - ResidualInftyErrPrev);
    fprintf('Iter %2d | Residual %2.5e \n',iter,ResiduaLInftyErr)
    if iter > 1e2
        break
    end
end

%%

LInftyErr1=max(abs(errUpper));
LInftyErr2=max(abs(errLower));
LInftyErr=[LInftyErr1, LInftyErr2];

L2Err1 = getL2normS4(dim, elemNonOverlap, center, h, n, errUpper);
L2Err2 = getL2normS4(dim, elemNonOverlap, center, h, n, errLower);
L2Err = [L2Err1, L2Err2];

EnergyErr1= sqrt(errUpper' * A * errUpper);
EnergyErr2= sqrt(errLower' * A * errLower);
EnergyErr=[EnergyErr1, EnergyErr2];


fprintf('L-infinity Error: [%e, %e]\n', LInftyErr1, LInftyErr2);
fprintf('L2 Error:         [%e, %e]\n', L2Err1, L2Err2);
fprintf('Energy Error:     [%e, %e]\n', EnergyErr1, EnergyErr2);