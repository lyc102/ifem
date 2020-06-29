%% SQUAREPOISSONVEM Poisson equation in a square domain.
%
%   squarePoissonVEM computes linear VEM approximations of the Poisson
%   equation in the unit square on a sequence of polygonal meshes generated
%   by PolyMesher.
% 
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
%
% Created by Min Wen.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


%% Load the information of mesh and PDE information
pde = vemcosdata;
%pde = vemdata;
nameV = [64,256,1024,4096, 20014]';
% nameV = [64,256,1024,4096]';
maxIt = length(nameV);
errorH1 = zeros(maxIt,1); %initialize the error
errorMax = zeros(maxIt,1); %initialize the error
solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1);
h = zeros(length(nameV),1); %initialize the error

%% Iterations
for k = 1:length(nameV)
    % load mesh
    load(['E',num2str(nameV(k)),'.mat']);
    % solve
    [u,A,asst,solt] = PoissonVEM(Node,Element,pde);    
    % plot
    if k <= 3 % only plot meshes of small size
        subplot(1,2,1); showmeshpoly(Node,Element);
        subplot(1,2,2); showsolutionpoly(Node,Element,u);
        pause(1);
    end
    assembleTime(k) = asst;
    solverTime(k) = solt;
    uI = pde.u(Node);
    errorH1(k) = sqrt((u-uI)'*A*(u-uI)); % get the error in H1 norm
    errorMax(k) = max(abs(u-uI));
    h(k) = 1/sqrt(length(u));
end

%% Plot convergence rates
figure;
showrateh2(h(1:k),errorH1(1:k),2,'r-+','||u_I-u_h||_A',...
           h(1:k),errorMax(1:k),2,'b-+','||u_I-u_h||_{\infty}');


%% Display error and time
err = struct('h',h,'elem',nameV,'H1',errorH1(1:k),'uIuhMax',errorMax(1:k));
time = struct('h',h,'solver',solverTime(1:k),'assemble',assembleTime(1:k));
disp('Table: CPU time')
colname = {'h','Assemble','Solve'};
displaytable(colname,time.h,'%0.3e',time.assemble,'%0.2e',time.solver,'%0.2e');
disp('Table: Error')
colname = {'h','NT','||DuI-Du_h||','||uI-u_h||_{max}'};     
displaytable(colname,err.h,'%0.3e',err.elem,[],err.H1,'%0.5e',err.uIuhMax,'%0.5e');
