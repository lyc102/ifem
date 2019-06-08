close all; clear;

%% Setting of the problem
global s
global mA

% problem set up 
s = 0.15;
h = 1/8;
maxIt = 3;
pde = fracLapdata2;
cube = [0 1 0 1 0 1];

% options for MG
mA = cell(maxIt+1,1);
option.tol = 1e-7;
option.gradmesh = 1;
option.smootherType = 0;
option.plotflag = 0;
option.smoothingstep = 3;

%% Finite element approximation
N = zeros(maxIt,1);
errH1 = zeros(maxIt,1);
for k = 1:maxIt
    tic;
    [uh,erru,errMG,node,elem] = mgfracLap2d(pde,cube,h,s,option);
    N(k) = length(uh);
    errH1(k) = erru;
    fprintf('#Dof %d   ',N(k));
    fprintf('MG converges at %d with error %e \n',length(errMG),errMG(end));
    toc;
    h = h/2;
end

%% Plot convergence rates
figure;
showrate(N,errH1,2,'-*');
display(errH1);