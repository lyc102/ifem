function [err,time,solver,eqn] = afemfracLap(node,elem,pde,bdFlag,option,varargin)
%% AFEMFRACLAP solve fractional Laplacian equation by various finite element methods
%
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if nargin >=1 && ischar(node)
    option.elemType = node;
    clear node
end
if ~exist('node','var') || ~exist('elem','var')
    [node,elem] = squaremesh([0,1,0,1],0.125);  % default mesh is a square
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = sincosdata;                          % default data
end
if ~exist('bdFlag','var')
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end

%% Parameters
option = afemoption(option,2);
if strcmp(option.elemType,'P1')
    option.elemType = 'P1P1';
end
maxIt = option.maxIt;   maxN = option.maxN; L0 = option.L0;
refType = option.refType; elemType = option.elemType; theta = option.theta;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
energy = zeros(maxIt,1);  errH1 = zeros(maxIt,1); 
errMax = zeros(maxIt,1);  erreta = zeros(maxIt,1); errestar = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
estimateTime = zeros(maxIt,1);
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    switch elemType
        case 'P1P1'     % P1 in x, P1 in y
            [u,eqn,info] = fracLapP1P1(node,elem,pde,option);
        case 'P1P2'     % P1 in x, P2 in y
            [u,eqn,info] = fracLapP1P2(node,elem,pde,option);
        case 'P2P2'     % P2 in x, P2 in y
            [u,eqn,info] = fracLapP2P2(node,elem,pde,option);
    end
    % compute error
    tic;
    if isfield(pde,'exactu')
        errH1(k) = info.errH1;
        % interpolation
        if strcmp(elemType(1:2),'P1')
            uI = Lagrangeinterpolate(pde.exactu,node,elem);
        elseif strcmp(elemType(1:2),'P2')
            uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
        end
        errMax(k) = max(abs(u(1:length(uI))-uI));
    end
    Nv = size(node,1);
    energy(k) = (u'*eqn.A*u)/2 - sum(eqn.b.*u);
    errTime(k) = toc;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),errH1(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(u);
    if option.plotflag && Nv < 3e3 % show mesh and solution for small size
       figure(1);  
       showresult(node,elem,u(1:Nv));  
    end
    if N(k) > maxN
        break;
    end
    % Step 2: ESTIMATE
    tic;
    switch option.estType
        case 'star'
            if isfield(option,'gNquadorder')
                option.gNquadorder = option.gNquadorder + 1;
            end
            [eta,estar] = estimatefracLap(node,elem,eqn.y,u,pde,option);
%             [eta,estar] = estimatefracLaposc(node,elem,eqn.y,u,pde,option);
        case 'recovery' % recovery type
            eta = estimaterecovery(node,elem,u(1:Nv));         
    end
    erreta(k) = sqrt(sum(eta.^2));
    errestar(k) = sqrt(sum(estar.^2)/3);
    estimateTime(k) = toc;
    % Step 3: MARK
    switch option.markType
        case 'L2'
            markedElem = mark(elem,eta,theta);
        case 'MAX'
            markedElem = mark(elem,eta,theta,'MAX');            
    end
    % Step 4: REFINE
    tic;
    if N(k) > option.maxN
        break;
    end
    [node,elem,bdFlag,HB,tree] = bisect(node,elem,markedElem,bdFlag); %#ok<ASGLU>
    % Step 4.2: COARSEN
    if option.coarsenflag 
        eta = eleminterpolate(eta,tree);
        markedElem = mark(elem,eta,0.25*theta,'COARSEN');
        [node,elem,bdFlag] = coarsen(node,elem,markedElem,bdFlag);
    end    
    meshTime(k) = toc;
end
% compute a fine energy
if ~isfield(pde,'exactu')
    [node,elem] = uniformrefine(node,elem);    
    [uh,eqn] = fracLapP1P1(node,elem,pde,option);
    fineEnergy = (uh'*eqn.A*uh)/2 - sum(eqn.b.*uh);
else
    fineEnergy = energy(k);
end
energyError = sqrt(2*(energy(1:k)-fineEnergy));

%% Plot convergence rates
if option.rateflag
    figure;
    if isfield(pde,'exactu')
        showrate3(N(1:k),errH1(1:k),k-9,'-*','||u-u_h||_A',...
                  N(1:k-1),erreta(1:k-1),k-9,'k-*','\eta',...
                  N(1:k-1),errestar(1:k-1),k-9,'r-*','E_star');
    else
        showrate3(N(1:k-2),energyError(1:k-2),8,'-*','||uf-uh||_A',...
                  N(1:k-1),erreta(1:k-1),k-9,'k-*','\eta',...
                  N(1:k-1),errestar(1:k-1),k-9,'r-*','E_star');
    end
end

%% Output
err = struct('N',N(1:k),'H1',errH1(1:k),'uIuhMax',errMax(1:k),...
             'energyError',energyError(1:k),'eta',erreta(1:k),'estar',errestar(1:k));
time = struct('N',N(1:k),'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k),...
              'estimate',estimateTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and time
fprintf('\n')
display('Table: Error')
if isfield(pde,'exactu')
    colname = {'#Dof','||u-uh||_A',' eta', ' Estar', '||uI-uh||_{max}','||uf-uh||_A'};
    disptable(colname,err.N,[],err.H1,'%0.5e',err.eta,'%0.5e',...
              err.estar,'%0.5e', err.uIuhMax,'%0.5e',err.energyError,'%0.5e');
else
    colname = {'#Dof','||uf-uh||_A',' eta', '  Estar'};
    disptable(colname,err.N,[],err.energyError,'%0.5e',err.eta,'%0.5e',err.estar,'%0.5e')
end

display('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Estimate','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.estimate,'%0.2e',time.mesh,'%0.2e');
% display('Aspect Ratio')
% fprintf('%4.2e \n',info.aspectRatio);
% figure; 
% plot(info.aspectRatio);
% hold on
% lr = length(info.aspectRatio);
% plot(1:lr,mean(info.aspectRatio)*ones(1,lr),'r')
display(mean(info.aspectRatio));              