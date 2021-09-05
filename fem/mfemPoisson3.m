function [err,time,solver,eqn] = mfemPoisson3(mesh,pde,option,varargin)
%% MFEMPOISSON3 solve Poisson equation by various mixed finite element methods
%
%   MFEMPOISSON computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
%
% See also femPoisson, mfemPoisson, Poisson3RT0, Poisson3BDM1
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = double(mesh.elem);
else
    [node,elem] = cubemesh([-1,1,-1,1,-1,1],1); % default mesh is a cube
end
if isfield(mesh,'bdFlag')
    bdFlag = mesh.bdFlag;
else
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = mixBCdata3;                          % default data
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
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
    end
end

%% Initialization
% error
erruL2 = zeros(maxIt,1);   
erruIuhL2 = zeros(maxIt,1); 
errsigmaL2 = zeros(maxIt,1); 
errsigmaHdiv = zeros(maxIt,1); 
errsigmaIsigmah = zeros(maxIt,1); 
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
        case 'RT0'  % RT0-P0 mixed FEM 
%             [u,sigma,eqn,info] = Poisson3RT0(node,elem,bdFlag,pde,option);
            [u,sigma,eqn,info] = Poisson3RT0(node,elem,bdFlag,pde,option);
        case 'BDM1' % BDM1-P0 mixed FEM
            [u,sigma,eqn,info] = Poisson3BDM1(node,elem,bdFlag,pde,option);
    end
    % compute error
    t = cputime;
    volume = simplexvolume(node,elem);
    if isfield(pde,'Du') && isfield(pde,'f')
        switch elemType
            case 'RT0'  % RT0-P0 mixed FEM 
                errsigmaL2(k) = getL2error3RT0(node,elem,pde.Du,sigma);
                divsigma = (eqn.B*sigma)./volume;
                errsigmaHdiv(k) = getL2error3(node,elem,pde.f,-divsigma);
                sigmaI = faceinterpolate3(pde.Du,node,elem);
            case 'BDM1' % BDM1-P0 mixed FEM
                errsigmaL2(k) = getL2error3BDM1(node,elem,pde.Du,sigma);
                divsigma = (eqn.B*sigma)./volume;
                errsigmaHdiv(k) = getL2error3(node,elem,pde.f,-divsigma);
                sigmaI = faceinterpolate3(pde.Du,node,elem,'BDM1');
        end
        errsigmaIsigmah(k)=sqrt((sigma-sigmaI)'*eqn.M*(sigma-sigmaI));
    end
    if isfield(pde,'exactu')
        erruL2(k) = getL2error3(node,elem,pde.exactu,u);
        % interpolation
        uI = Lagrangeinterpolate(pde.exactu,node,elem,'P0');
        erruIuhL2(k) = sqrt(dot((uI-u).^2,volume));
    end
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
            errTime(k),erruL2(k), errsigmaL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = size(node,1);
    h(k) = 1./(sqrt(size(node,1))-1);    
    if option.plotflag && N(k) < 2e3 % show mesh and solution for small size
       figure(1);  
%        showresult3(node,elem,u);    
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
    showrateh2(h(1:k),erruIuhL2(1:k),1,'-*','|| u_I-u_h||_{\infty}',...
               h(1:k),erruL2(1:k),1,'k-+','|| u-u_h||');
    subplot(1,2,2)
    showrateh3(h(1:k),errsigmaHdiv(1:k),1,'-+','|| div(\sigma - \sigma_h)||',... 
               h(1:k),errsigmaL2(1:k),1,'k-*','|| \sigma - \sigma_h||',...
               h(1:k),errsigmaIsigmah(1:k),1,'m-+','|| \sigma_I - \sigma_h||');
end

%% Output
err = struct('h',h(1:k),'N',N(1:k),'uL2',erruL2(1:k),'uIuhL2',erruIuhL2(1:k),...
             'sigmaL2',errsigmaL2(1:k),'sigmaHdiv',errsigmaHdiv(1:k),...
             'sigmaIsigmahL2',errsigmaIsigmah(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and time
if option.dispflag
    fprintf('\n');
%     disp('Table: Error')
    colname = {'#Dof','h','||u-u_h||','||u_I-u_h||','||sigma-sigma_h||','||sigma-sigma_h||_{div}'};
    displaytable(colname,err.N,[],err.h,'%0.2e',err.uL2,'%0.5e',err.uIuhL2,'%0.5e',...
                     err.sigmaL2,'%0.5e',err.sigmaHdiv,'%0.5e');
                 
%     disp('Table: CPU time')
    colname = {'#Dof','Assemble','Solve','Error','Mesh'};
    disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');   
    fprintf('\n');
end