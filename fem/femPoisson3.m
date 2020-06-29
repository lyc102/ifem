function [err,time,solver,eqn] = femPoisson3(mesh,pde,option,varargin)
%% FEMPOISSON3 solve Poisson equation by various finite element methods
%
%   FEMPOISSON3 computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% See also Poisson, crack, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = double(mesh.elem);
else
    [node,elem] = cubemesh([0,1,0,1,0,1],0.25); % default mesh is a cube
end
if isfield(mesh,'bdFlag')
    bdFlag = mesh.bdFlag;
else
    bdFlag = setboundary3(node,elem,'Dirichlet'); 
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = sincosdata3;                          % default data
end

%% Parameters
option = femoption(option);
maxIt = option.maxIt;   
maxN = option.maxN; 
L0 = option.L0;
elemType = option.elemType; 
refType = option.refType;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
end

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1); h = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
%     bdFlag = sparse(double(bdFlag));
    % solve the equation
    switch elemType
        case 'P1'     % piecewise linear function P1 element
            [soln,eqn,info] = Poisson3(node,elem,bdFlag,pde,option);
        case 'CR'     % piecewise linear function CR element
            [soln,eqn,info] = Poisson3CR(node,elem,bdFlag,pde,option);
        case 'P2'     % piecewise quadratic function
            [soln,eqn,info] = Poisson3P2(node,elem,bdFlag,pde,option);
        case 'WG'     % weak Galerkin element
            [soln,eqn,info] = Poisson3WG(node,elem,bdFlag,pde,option);            
    end
    uh = soln.u;
    % compute error
    t = cputime;
    if isfield(pde,'Du')
        if isfield(soln,'Du') && ~isempty(soln.Du) % Du is in the output
            errH1(k) = getH1error3(node,elem,pde.Du,soln.Du);
        else
            errH1(k) = getH1error3(node,elem,pde.Du,uh);            
        end
    end
    if isfield(pde,'exactu')
        errL2(k) = getL2error3(node,elem,pde.exactu,uh);        
        % interpolation
        switch elemType
            case 'P1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'CR'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',eqn.face);
            case 'P2'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
            case 'WG'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'WG',eqn.face);
        end
        erruIuh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
        errMax(k) = max(abs(uh-uI));
    end
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
                 errTime(k), errH1(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(soln.u);
    h(k) = 1./(size(node,1)^(1/3)-1);    
    if  strcmp(elemType,'WG') % modify size for WG
        if ~isfield(option,'reducesystem') || (option.reducesystem == 1)
            N(k) = N(k) - size(elem,1); % reduced system
        end    
    end                
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
        switch elemType
        case 'P1'     % piecewise linear function P1 element
            figure(1);  showresult3(node,elem,uh);    
        case 'CR'     % piecewise linear function CR element
            continue;
        case 'P2'     % piecewise quadratic function
            figure(1);  showresult3(node,elem,uh(1:size(node,1)));    
        case 'WG'     % weak Galerkin element
            continue;
        end
    end
    if N(k) > maxN
        break;
    end
    % refine mesh
    t = cputime;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
    meshTime(k) = cputime - t;
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh2(h(1:k),errH1(1:k),1,'-*','||Du-Du_h||',...
               h(1:k),errL2(1:k),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrateh2(h(1:k),erruIuh(1:k),1,'m-+','||Du_I-Du_h||',...
               h(1:k),errMax(1:k),1,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('h',h(1:k),'N',N,'H1',errH1(1:k),'L2',errL2(1:k),...
             'uIuhH1',erruIuh(1:k),'uIuhMax',errMax(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));

%% Display error and time
if option.dispflag
    disp('Table: Error')
    colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||DuI-Du_h||','||uI-u_h||_{max}'};
    disptable(colname,err.N,[],err.h,'%0.3e',err.L2,'%0.5e',err.H1,'%0.5e',...
              err.uIuhH1,'%0.5e',err.uIuhMax,'%0.5e');

    disp('Table: CPU time')
    colname = {'#Dof','Assemble','Solve','Error','Mesh'};
    disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                      time.err,'%0.2e',time.mesh,'%0.2e');
end
