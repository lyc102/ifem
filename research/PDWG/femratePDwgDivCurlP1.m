%% PDWG rate for P1
clear; close all;
%%
maxIt = 4;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2UTotal = zeros(maxIt,1);
errqlambdaTotal = zeros(maxIt,1);
errSTotal = zeros(maxIt,1);
errL2QuTotal = zeros(maxIt,1);
%% data
pde = DataDivCurlSingular3;

%% Finite Element Method     
for k = 1:maxIt
    h(k) = (0.5)^k;
    [node,elem] = cubemesh([0,1,0,1,0,1],h(k));
    bdFlag = setboundary3(node,elem,'Dirichlet');
    [soln,eqn,info] = DivCurl3PDWGP1(node,elem,bdFlag,pde);
    N(k) = size(elem,1);
    
    errL2UTotal(k) = eqn.errorU;
    errqlambdaTotal(k) = eqn.errorqlambda;
    errSTotal(k) = eqn.errorS;
    [errL2QuTotal(k), ~] = getL2error3ND1d(node,elem,soln.u,pde);
    fprintf('PDWG method with h      %3g \n', h(k));
    fprintf('Number of Dof           %d \n', sum(eqn.freeDof));
    fprintf('L2 error in u:          %6g \n', errL2UTotal(k));
    fprintf('L2 error in (q,lambda): %6g \n', errqlambdaTotal(k));
    fprintf('L2 error in s:          %6g \n', errSTotal(k));
    fprintf('L2 error in Qu-u_h:     %6g \n\n', errL2QuTotal(k));   
end

%%
visualizationPDwgDivCurl;