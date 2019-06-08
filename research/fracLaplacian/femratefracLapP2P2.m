close all; clear all;

%% Setting of the problem
global s
s = 0.7;
maxIt = 4;
pde = fracLapdata2;
pde.s = s;
pde.L = 1;
option.solver = 'direct';
% option.solver = 'amg';
option.gNquadorder = 4;
option.tol = 1e-6;
[node,elem] = squaremesh([0 1 0 1],1);

%% Finite element approximation
N = zeros(maxIt,1);
errH1 = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);
for k = 1:maxIt
    [node,elem] = uniformrefine(node,elem);
    tic;
    [uh,eqn,info] = fracLapP2P2(node,elem,pde,option);
    N(k) = length(uh);
    uI = pde.exactu(node); % nodal interpolation
%     figure(1); showsolution(node,elem,uI);
%     figure(2); showsolution(node,elem,uh(1:size(node,1)));
%     erruIuh(k) = norm(uI-uh(1:size(node,1)))/size(elem,1);
    erruIuh(k) = max(abs(uI-uh(1:size(node,1))));
%     erruIuh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
%     errH1(k) = getH1error3bd(node,elem,pde,uh,eqn.A,5);
    errH1(k) = info.errH1;
    fprintf('#Dof %d   %0.5e',N(k),errH1(k));    
    toc;
end

%% Plot convergence rates
figure(2)
showrate2(N,errH1,2,'-*','||u-u_h||_A',N,erruIuh,2,'r-*','||u_I-u_h||_{\infty}');

%% Display error
colname = {' #Dof   ||Du-Du_h||_w   ||u_I-u_h||_{\infty}'};
disptable(colname,N,[],errH1,'%0.5e',erruIuh,'%0.5e')