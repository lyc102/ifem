function [err,time,solver,eqn] = femfracLap(node,elem,pde,bdFlag,option,varargin)
%% FEMFRACLAP solve fractional Laplacian equation by various finite element methods
%
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

nv = size(elem,2);
% dim = size(node,2);

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
option = femoption(option);
if strcmp(option.elemType,'P1')
    option.elemType = 'P1P1';
end
maxIt = option.maxIt;   maxN = option.maxN; L0 = option.L0;
refType = option.refType; elemType = option.elemType; 

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        if nv == 4
            [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
        else
            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        end        
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
energy = zeros(maxIt,1);  errH1 = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch elemType
        case 'P1P1'     % P1 in x, P1 in y
            [u,eqn,info] = fracLapP1P1(node,elem,pde,option);
        case 'P1P2'     % P1 in x, P2 in y
            [u,eqn,info] = fracLapP1P2(node,elem,pde,option);
        case 'P2P2'     % P2 in x, P2 in y
            [u,eqn,info] = fracLapP2P2(node,elem,pde,option);
        case 'P2bP2'     % P2b in x, P2 in y
            [u,eqn,info] = fracLapP2bP2(node,elem,pde,option);
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
       figure(1);  showresult(node,elem,u(1:Nv));  
    end
    if N(k) > maxN
        break;
    end
    % refine mesh
    tic;
    if strcmp(refType,'red')
        if nv == 4
            [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
        else
            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        end        
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    meshTime(k) = toc;
end
% compute a fine energy
if ~isfield(pde,'exactu')
%     [node,elem] = uniformrefine(node,elem);    
    switch elemType
        case 'P1P1'     % P1 in x, P1 in y
            [u,eqn] = fracLapP1P1(node,elem,pde,option);
        case 'P1P2'     % P1 in x, P2 in y
            [u,eqn] = fracLapP1P2(node,elem,pde,option);
        case 'P2P2'     % P2 in x, P2 in y
            [u,eqn] = fracLapP2P2(node,elem,pde,option);
        case 'P2bP2'     % P2b in x, P2 in y
            [u,eqn] = fracLapP2bP2(node,elem,pde,option);            
    end
    fineEnergy = (u'*eqn.A*u)/2 - sum(eqn.b.*u);
else
    fineEnergy = energy(k);
end
energyError = sqrt(2*(energy(1:k)-fineEnergy));

%% Plot convergence rates
if option.rateflag
    figure;
    if isfield(pde,'exactu')
        showrate(N(1:k),errH1(1:k),2,'-*','||u-u_h||_A');
    else
        showrate(N(1:k),energyError(1:k),2,'r-*','||uf-uh||_A');
    end
end

%% Output
err = struct('N',N(1:k),'H1',errH1(1:k),'uIuhMax',errMax(1:k),...
             'energyError',energyError(1:k));
time = struct('N',N(1:k),'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and time
fprintf('\n')
display('Table: Error')
if isfield(pde,'exactu')
    colname = {'#Dof','||u-u_h||_A','||uI-u_h||_{max}','||uf-uh||_A'};
    disptable(colname,err.N,[],err.H1,'%0.5e',err.uIuhMax,'%0.5e',energyError,'%0.5e');
else
    colname = {' #Dof','||uf-uh||_A'};
    disptable(colname,N,[],energyError,'%0.5e')
end

display('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');