function [err,time,solver,eqn] = femHelmholtz(node,elem,pde,bdFlag,option,varargin)
%% FEMPOISSON solve Poisson equation by various finite element methods
%
%   FEMPOISSON computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% See also Poisson, crack, Lshape
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
option = femoption(option);
maxIt = option.maxIt;   maxN = option.maxN; L0 = option.L0;
elemType = option.elemType; refType = option.refType;

%% Generate an initial mesh 
for i = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Finite Element Method        
for i = 1:maxIt
    % solve the equation
    switch elemType
        case 'P1'     % piecewise linear function P1 element
            [u,eqn,info] = Helmholtz(node,elem,pde,bdFlag,option);
            u = real(u);
        case 'P2'     % piecewise quadratic function
            [u,eqn,info] = HelmholtzP2(node,elem,pde,bdFlag,option);
        case 'P3'     % piecewise cubic function
            [u,eqn,info] = HelmholtzP3(node,elem,pde,bdFlag,option);
        case 'WG'     % weak Galerkin element
            [u,eqn,info] = HelmholtzWG(node,elem,pde,bdFlag,option);            
    end
    % compute error
    tic;
    if isfield(pde,'Du') 
        if exist('Du','var') && ~isempty(Du) 
            errH1(i) = getH1error(node,elem,pde.Du,Du);
        else
            errH1(i) = getH1error(node,elem,pde.Du,u);            
        end
    end
    if isfield(pde,'exactu')
        errL2(i) = getL2error(node,elem,pde.exactu,u);
        % interpolation
        if strcmp(elemType,'P1')
            uI = Lagrangeinterpolate(pde.exactu,node,elem);
        else
            uI = Lagrangeinterpolate(pde.exactu,node,elem,elemType,eqn.edge);
        end
        erruIuh(i) = sqrt((u-uI)'*(eqn.Delta)*(u-uI));
        errMax(i) = max(abs(u-uI));
    end
    errTime(i) = toc;
    % record time
    solverTime(i) = info.solverTime;
    assembleTime(i) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(i),errH1(i), errL2(i));    
    end
    % record solver information
    itStep(i) = info.itStep;
    stopErr(i) = info.stopErr;
    flag(i) = info.flag;
    % plot 
    N(i) = length(u);
    if option.plotflag && N(i) < 2e3 % show mesh and solution for small size
       figure(1);  showresult(node,elem,u);    
    end
    if N(i) > maxN
        break;
    end
    % refine mesh
    tic;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    meshTime(i) = toc;
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrate2(N(1:i),errH1(1:i),1,'-*','||Du-Du_h||',...
              N(1:i),errL2(1:i),1,'k-+','||u-u_h||');
    subplot(1,2,2)
    showrate2(N(1:i),erruIuh(1:i),1,'m-+','||Du_I-Du_h||',...
              N(1:i),errMax(1:i),1,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('N',N(1:i),'H1',errH1(1:i),'L2',errL2(1:i),...
             'uIuhH1',erruIuh(1:i),'uIuhMax',errMax(1:i));
time = struct('N',N(1:i),'err',errTime(1:i),'solver',solverTime(1:i), ...
              'assmble',assembleTime(1:i),'mesh',meshTime(1:i));
solver = struct('N',N(1:i),'itStep',itStep(1:i),'time',solverTime(1:i),...
                'stopErr',stopErr(1:i),'flag',flag(1:i));
            
%% Display error and iteration 
global k     
n = floor(sqrt(err.N)-1);
colname = {'#n','||u-u_h||','iter', 'solvetime','k*h'};
display(colname, n,[],err.L2,'%0.5e',solver.itStep,[],...
                 solver.time,'%0.2e',k./n,'%0.3e');     