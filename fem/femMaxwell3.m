function [err,time,solver,eqn] = femMaxwell3(mesh,pde,option,varargin)
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
    pde = Maxwelldata2;                          % default data
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
errL2 = zeros(maxIt,1);   errHcurl = zeros(maxIt,1); 
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
        case 'ND0'     % the lowest order edge element
            if(abs(pde.epsilon)>1.0e-8)
                [u,edge,eqn,info] = Maxwell(node,elem,bdFlag,pde,option);
            else
                [u,edge,eqn] = Maxwellsaddle(node,elem,bdFlag,pde,option); 
            end                
        case 'ND1'     % the lowest order second family
            if(abs(pde.epsilon)>1.0e-8)
                [u,edge,eqn,info] = Maxwell1(node,elem,bdFlag,pde,option);
            else
                [u,edge,eqn] = Maxwellsaddle1(node,elem,bdFlag,pde,option); 
            end                
        case 'ND2'     % quadratic Nedelec element
            if(abs(pde.epsilon)>1.0e-8)
                [u,T,eqn,info] = Maxwell2(node,elem,bdFlag,pde,option); 
            else
                [u,edge,eqn] = Maxwellsaddle2(node,elem,bdFlag,pde,option); 
            end                
    end
    % compute error
    t = cputime;
    if isfield(pde,'curlu')
        switch elemType
            case 'ND0'            
               errHcurl(k) = getHcurlerror3ND(node,elem,pde.curlu,real(u));
            case 'ND1'
               errHcurl(k) = getHcurlerror3ND1(node,elem,pde.curlu,real(u));
            case 'ND2'
               errHcurl(k) = getHcurlerror3ND2(node,elem,pde.curlu,real(u));               
        end
    end
    if isfield(pde,'exactu')
        % interpolation
        switch elemType
            case 'ND0'
                errL2(k) = getL2error3ND(node,elem,pde.exactu,real(u));        
                uI = edgeinterpolate(pde.exactu,node,edge);
            case 'ND1'
                errL2(k) = getL2error3ND1(node,elem,pde.exactu,real(u));        
                uI = edgeinterpolate1(pde.exactu,node,edge);
            case 'ND2'
                errL2(k) = getL2error3ND2(node,elem,pde.exactu,real(u));  
                uI = edgeinterpolate2(pde.exactu,node,T.edge,T.face,T.face2edge);
        end
        erruIuh(k) = sqrt((u-uI)'*eqn.A*(u-uI));
        errMax(k) = max(abs(u-uI));
    end
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
                 errTime(k), errHcurl(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);    
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
    showrateh2(h(1:k),errHcurl(1:k),1,'-*','||Du-Du_h||',...
               h(1:k),errL2(1:k),1,'k-+','||u-u_h||');
    title(['Convergence Rate of 3D ', elemType, '-element']);
    subplot(1,2,2)
    showrateh2(h(1:k),erruIuh(1:k),1,'m-+','||Du_I-Du_h||',...
               h(1:k),errMax(1:k),1,'r-*','||u_I-u_h||_{\infty}');
    title(['Convergence Rate of 3D ', elemType, '-element']);
end

%% Output
err = struct('h',h(1:k),'N',N,'H1',errHcurl(1:k),'L2',errL2(1:k),...
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
