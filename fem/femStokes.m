function [err,time,solver,eqn] = femStokes(mesh,pde,option,varargin)
%% FEMSTOKES solve the Stokes equation by various finite element methods
%
%   FEMSTOKES computes approximations to the Stokes equation on a
%   sequence of meshes obtained by uniform refinement of the input mesh.
% 
% See also femPoisson, femStokesHdiv
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = double(mesh.elem);
else
    [node,elem] = squaremesh([0,1,0,1],0.125);  % default mesh is a square
end
if isfield(mesh,'bdFlag')
    bdFlag = mesh.bdFlag;
else
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = Stokesdata1;                          % default data
end

%% Parameters
option = femoption(option);
maxIt = option.maxIt;   
maxN = option.maxN; 
L0 = option.L0;
option = femStokesoption(option);
elemType = option.elemType; 
refType = option.refType;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(option.refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(option.refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
erruH1 = zeros(maxIt,1); errpL2 = zeros(maxIt,1); 
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
    t = cputime;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    meshTime(k) = cputime - t;    
    % solve the equation
    switch upper(elemType)
        case 'P2P1'
            [soln,eqn,info] = StokesP2P1(node,elem,bdFlag,pde,option);
        case 'P2P0'
            [soln,eqn,info] = StokesP2P0(node,elem,bdFlag,pde,option);
        case 'ISOP2P1'
            [soln,eqn,info] = StokesisoP2P1(node,elem,bdFlag,pde,option);
        case 'ISOP2P0'
            [soln,eqn,info] = StokesisoP2P0(node,elem,bdFlag,pde,option);
        case 'CRP0'
            [soln,eqn,info] = StokesCRP0(node,elem,bdFlag,pde,option);
        case 'CRP1'
            [soln,eqn,info] = StokesCRP1(node,elem,bdFlag,pde,option);
%         case 'MINI'
%             [soln,eqn,info] = StokesMini(node,elem,bdFlag,pde,option);
        case 'P1BP1'
            [soln,eqn,info] = StokesP1bP1(node,elem,bdFlag,pde,option);
        case 'RTP0'
            [soln,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option);        
        case 'BDMP0'
            [soln,eqn,info] = StokesBDM1B(node,elem,bdFlag,pde,option);                    
    end    
    % compute error
    t = cputime;
    uh = soln.u;
    ph = soln.p;
    if strcmp(elemType,'P2P0') || strcmp(elemType,'P2P1') || ...
       strcmp(elemType,'isoP2P0') || strcmp(elemType,'isoP2P1')
        uI = pde.exactu([node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2]);
    elseif strcmp(elemType,'CRP0') || strcmp(elemType,'CRP1')
        uI = pde.exactu((node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2);
%     elseif strcmp(elemType,'Mini')
%         uI = pde.exactu(node);
    elseif strcmp(elemType,'P1bP1')
        % bubble part won't be taken in the error computation
        Nv = size(node,1); NT = size(elem,1);
        u0 = pde.exactu(node);
        uI = uh;
        uI([(1:Nv)'; NT+Nv+(1:Nv)']) = u0(:);
    end
    erruH1(k) = sqrt((uh-uI(:))'*eqn.Lap*(uh-uI(:)));
    errpL2(k) = getL2error(node,elem,pde.exactp,ph);
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),erruH1(k), errpL2(k));
    end
    % record solver information
    itStep(k) = info.itStep;
%    stopErr(k) = info.stopErr;
%    flag(k) = info.flag;
    % plot
    N(k) = size(node,1);
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
        if length(ph) == size(elem,1) % piecewise constant function
            p1 = recoverP02P1(node,elem,ph); % recovery a linear velocity
        elseif length(ph) == size(node,1)
           p1 = ph;
        end
        figure(1);  showresult(node,elem,p1);
    end
    if N(k) > maxN
        break;
    end
end
h = 1./(sqrt(N(1:k))-1);

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh(h,erruH1(1:k),1,'k-+','|u_I-u_h|_1');
    subplot(1,2,2)
    showrateh(h,errpL2(1:k),1,'m-+','|| p-p_h||');
end

%% Output
err = struct('h',h,'N',N,'uH1',erruH1(1:k),'pL2',errpL2(1:k));          
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error
disp('Table: Error')
colname = {'#nodes'     '|u_I-u_h|_1'      '||p-p_h||'};
disptable(colname,err.N,[],err.uH1,'%0.5e',err.pL2,'%0.5e');

disp('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');
