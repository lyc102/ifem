close all; clear all;

%% Setting of the problem
global s
s = 0.3;
maxIt = 4;
pde = fonedata;
pde.s = s;
pde.L = 1;
option.solver = 'direct';
% option.solver = 'amg';
option.gNquadorder = 4;
[node,elem] = squaremesh([-1,1,-1,1],2);
% [node,elem] = delmesh(node,elem,'x>0 & y<0');

%% Finite element approximation
N = zeros(maxIt,1);
energy = zeros(maxIt,1);
for k = 1:maxIt
    [node,elem] = uniformrefine(node,elem);
    tic;
    [uh,eqn] = fracLapP1P1(node,elem,pde,option);
    Nv = size(node,1);
    figure(1); showresult(node,elem,uh(1:Nv)); 
    pause(0.1)     
    N(k) = length(uh);
    energy(k) = (uh'*eqn.A*uh)/2 - sum(eqn.b(1:Nv).*uh(1:Nv));
    fprintf('#Dof %d   ',N(k));    
    toc;
end
[node,elem] = uniformrefine(node,elem);    
tic;
[uh,eqn,info] = fracLapP1P1(node,elem,pde,option);
toc;
display(length(uh));
fineEnergy = (uh'*eqn.A*uh)/2 - sum(eqn.b.*uh);

%% Plot convergence rates
energyError = sqrt(2*(energy(1:k)-fineEnergy));
figure(2)
showrate(N(1:k),energyError(1:k),1,'r-*','Energy Error');

%% Display error
colname = {' #Dof   sqrt(2*(E(uh)-E(u)))'};
disptable(colname,N,[],energyError,'%0.5e')