close all; clear all;

%% Setting of the problem
global s
s = 0.3;
maxIt = 50;
maxN = 2e4;
theta = 0.3;
pde = fonedata;
pde.s = s;
pde.L = 1;
option.solver = 'direct';
% option.solver = 'amg';
option.gNquadorder = 4;
[node,elem] = squaremesh([-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x>0 & y<0');

%% Finite element approximation
N = zeros(maxIt,1);
energy = zeros(maxIt,1);
for k = 1:maxIt
    tic;
    % Step 1: SOLVE
    [uh,eqn] = fracLapP1P1(node,elem,pde,option);
    Nv = size(node,1);
    figure(1); showresult(node,elem,uh(1:Nv)); 
    pause(0.1) 
    N(k) = length(uh);
    energy(k) = (uh'*eqn.A*uh)/2 - sum(eqn.b(1:Nv).*uh(1:Nv));
    fprintf('#Dof %d   ',N(k));    
    toc;
    if N(k) > maxN
        break
    end
    % Step 2: ESTIMATE
%     eta = estimatefracLap(node,elem,eqn.y,uh,pde,option);
    eta = estimaterecovery(node,elem,uh(1:Nv));
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem] = bisect(node,elem,markedElem);    
end
figure(2); showmesh(node,elem);
[node,elem] = uniformrefine(node,elem);    
tic;
[uh,eqn,info] = fracLapP1P1(node,elem,pde,option);
toc;
display(length(uh));
fineEnergy = (uh'*eqn.A*uh)/2 - sum(eqn.b.*uh);

%% Plot convergence rates
energyError = sqrt(2*(energy(1:k)-fineEnergy));
figure(3)
showrate(N(1:k-1),energyError(1:k-1),5,'r-*','Energy Error');

%% Display error
colname = {' #Dof   sqrt(2*(E(uh)-E(u)))'};
disptable(colname,N(1:k),[],energyError(1:k),'%0.5e')