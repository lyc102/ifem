close all; clear all;

%% Setting of the problem
global s
s = 0.2;
maxIt = 4;
h = 1/4;
pde = fracLapdata2;
option.solver = 'direct';
option.gNquadorder = 4;
option.gradmesh = 1;
option.tol = 1e-6;

cube = [0,1,0,1,0,1];

%% Finite element approximation
N = zeros(maxIt,1);
errH1 = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);
for k = 1:maxIt
%     cube = [0,1,0,1,0,floor(abs(log(h)))];
    tic;
    [uh,eqn,info,node,elem,Neumann] = fracLap2d(cube,h,pde,option);
    N(k) = length(uh);
    fprintf('#Dof %d   ',N(k));    
    toc;
    uI = pde.exactu(node); % nodal interpolation
    erruIuh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
%     errH1(k) = getH1error3Q1(node,elem,pde.Du,uh,pde.d);       
    errH1(k) = info.errH1;
    h = h/2;
end

%% Plot convergence rates
figure(2)
showrate2(N,errH1,2,'-*','||u-u_h||_A',N,erruIuh,2,'r-*','||u_I-u_h||_A');

%% Display error
ts = zeros(k,2); ts = char(ts);
display(' #Dof   ||Du-Du_h||_w   ||DuI-Du_h||_w');
display([num2str(N) ts num2str(errH1,'%0.5e') ts num2str(erruIuh,'%0.5e')]);