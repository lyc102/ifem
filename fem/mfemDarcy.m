function [err,time,solver,eqn] = mfemDarcy(mesh,pde,option,varargin)
%% MFEMDARCY solve Poisson equation by various mixed finite element methods
%
%   MFEMDARCY computes approximations to the Darcy equation on a
%   sequence of meshes obtained by uniform refinement of an input mesh.
% 
%   It is almost identitical to mfemPoisson except the name of unknown is
%   changed to (u,p).
%
% See also mfemPoisson
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
    pde = mixBCdata;                          % default data
end

%% Parameters
option = mfemoption(option);
elemType = option.elemType;     
refType = option.refType;
maxIt = option.maxIt;       
maxN = option.maxN;         
L0 = option.L0;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    if isfield(pde,'K') && isnumeric(pde.K)
        pde.K = repmat(pde.K,4,1); % prolongate piecwise K to the fine grid
    end    
end

%% Initialize err and time array
% error
errpL2 = zeros(maxIt,1);   
errpIphL2 = zeros(maxIt,1); 
erruL2 = zeros(maxIt,1);   
erruHdiv = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); 
% time
errTime = zeros(maxIt,1);      
solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); 
meshTime = zeros(maxIt,1); 
% info
itStep = zeros(maxIt,1);  
stopErr = zeros(maxIt,1); 
flag = zeros(maxIt,1);
N = zeros(maxIt,1); 
h = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    switch elemType
        case 'RT0'  % RT0 mixed FEM 
            [p,u,eqn,info] = DarcyRT0(node,elem,bdFlag,pde,option);
        case 'BDM1' % BDM1 mixed FEM
            [p,u,eqn,info] = DarcyBDM1(node,elem,bdFlag,pde,option);
    end
    % compute error
    t = cputime;
    if isfield(pde,'exactu') && isfield(pde,'f')
        if strcmp(elemType,'RT0')
            erruL2(k) = getL2errorRT0(node,elem,pde.exactu,u);
            erruHdiv(k) = getHdiverrorRT0(node,elem,pde.f,-u,[]);
            uI = faceinterpolate(pde.exactu,node,elem,'RT0');
        else
            erruL2(k) = getL2errorBDM1(node,elem,pde.exactu,u);
            erruHdiv(k) = getHdiverrorBDM1(node,elem,pde.f,-u,[]);
            uI = faceinterpolate(pde.exactu,node,elem,'BDM1');
        end
        erruIuh(k)=sqrt((u-uI)'*eqn.M*(u-uI));
    end
    if isfield(pde,'exactp')
        errpL2(k) = getL2error(node,elem,pde.exactp,p);
        % interpolation
        pI = Lagrangeinterpolate(pde.exactp,node,elem,'P0');
        area = simplexvolume(node,elem);
        errpIphL2(k) = sqrt(dot((pI-p).^2,area));
    end
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),errpL2(k), erruL2(k));
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(p) + length(u);
    h(k) = 1./(sqrt(size(node,1))-1);
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
       figure(1);  
       showresult(node,elem,p);    
    end
    if N(k) > maxN
        break;
    end
    % refine mesh
    t = cputime;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    if isfield(pde,'K') && isnumeric(pde.K)
        pde.K = repmat(pde.K,4,1); % prolongate to the fine grid
    end
    meshTime(k) = cputime - t;
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh2(h(1:k),errpIphL2(1:k),1,'-*','||p_I-p_h||_{\infty}',...
               h(1:k),errpL2(1:k),1,'k-+','||p-p_h||');
    subplot(1,2,2)
    showrateh3(h(1:k),erruHdiv(1:k),1,'-+','||div(u - u_h)||',... 
               h(1:k),erruL2(1:k),1,'k-*','|| u - u_h||',...
               h(1:k),erruIuh(1:k),1,'m-+','|| u_I - u_h||');
end

%% Output
err = struct('h',h(1:k),'N',N(1:k),'pL2',errpL2(1:k),'pIphL2',errpIphL2(1:k),...
             'uL2',erruL2(1:k),'uHdiv',erruHdiv(1:k),...
             'uIuL2',erruIuh(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and CPU time
disp('Table: Error')
colname = {'#Dof','h','||p-p_h||','||p_I-p_h||','||u-u_h||','||u-u_h||_{div}'};
disptable(colname,err.N,[],err.h,'%0.2e',err.pL2,'%0.5e',err.pIphL2,'%0.5e',...
                     err.uL2,'%0.5e',err.uHdiv,'%0.5e');
                 
disp('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');     