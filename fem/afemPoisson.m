function [err,time,solver,eqn,node,elem] = afemPoisson(mesh,pde,option,varargin)

%% Check input arguments
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = double(mesh.elem);
end
if ~exist('node','var') || ~exist('elem','var')
    % default mesh: Lshape
    [node,elem] = squaremesh([-1,1,-1,1],1);
    [node,elem] = delmesh(node,elem,'x>0 & y<0');
end
if isfield(mesh,'bdFlag')
    bdFlag = mesh.bdFlag;
else % default boundary condition
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = Lshapedata;                          % default data
end

%% Parameters
option = afemoption(option,2);
maxIt = option.maxIt;
refType = option.refType;
elemType = option.elemType;
theta = option.theta;

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); erreta = zeros(maxIt,1);
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1);
estimateTime = zeros(maxIt,1);
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Generate an initial mesh 
for k = 1:option.L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    switch elemType
        case 'P1'     % piecewise linear function P1 element
            [soln,eqn,info] = Poisson(node,elem,bdFlag,pde,option);
        case 'CR'     % piecewise linear function CR element
            [soln,eqn,info] = PoissonCR(node,elem,bdFlag,pde,option);
        case 'P2'     % piecewise quadratic function
            [soln,eqn,info] = PoissonP2(node,elem,bdFlag,pde,option);
        case 'WG'     % weak Galerkin element
            [soln,eqn,info] = PoissonWG(node,elem,bdFlag,pde,option);            
    end
    uh = soln.u;
    % compute error
    t = cputime;
    if isfield(pde,'Du')
        if ~isfield(pde,'d')
            pde.d = [];
        end
        if isfield(soln,'Du') && ~isempty(soln.Du) % Du is in the output
            errH1(k) = getH1error(node,elem,pde.Du,soln.Du,pde.d);
        else
            errH1(k) = getH1error(node,elem,pde.Du,soln.Du,pde.d);            
        end
    end
    if isfield(pde,'exactu')        
        errL2(k) = getL2error(node,elem,pde.exactu,uh);
        % interpolation
        switch elemType
            case 'P1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'Q1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'CR'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',eqn.edge);
            case 'P2'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
            case 'P3'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P3',eqn.edge);
            case 'WG'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'WG',eqn.edge);
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
            errTime(k),errH1(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(uh);
    if  strcmp(elemType,'WG') % modify size for WG
        if ~isfield(option,'reducesystem') || (option.reducesystem == 1)
            N(k) = N(k) - size(elem,1); % reduced system
        end    
    end                
    if option.plotflag % show mesh and solution
       figure(1);  
       showresult(node,elem,uh,option.viewangle);   
    end
    % Step 2: ESTIMATE
    t = cputime;
    switch option.estType
        case 'recovery' % recovery type
            eta = estimaterecovery(node,elem,uh);         
        case 'residual' % residual type
            switch elemType
                case 'P1'
                    eta = estimateresidual(node,elem,uh,pde,bdFlag);
                case 'WG'
                    eta = estimateresidualWG(node,elem,uh,soln.Du,pde);                    
            end
    end
    erreta(k) = sqrt(sum(eta.^2));
    estimateTime(k) = cputime - t;
    % Step 3: MARK
    switch option.markType
        case 'L2'
            markedElem = mark(elem,eta,theta);
        case 'MAX'
            markedElem = mark(elem,eta,theta,'MAX');            
    end
    % Step 4: REFINE
    t = cputime;
    if (N(k) > option.maxN) || (k == maxIt) % no refinement
        break;
    end
    [node,elem,bdFlag,HB,tree] = bisect(node,elem,markedElem,bdFlag); %#ok<ASGLU>
    % Step 4.2: COARSEN
    if option.coarsenflag 
        eta = eleminterpolate(eta,tree);
        markedElem = mark(elem,eta,0.25*theta,'COARSEN');
        [node,elem,bdFlag] = coarsen(node,elem,markedElem,bdFlag);
    end
    meshTime(k) = cputime - t;
end

%% Plot convergence rates
if option.rateflag
    if ~isfield(option,'rateshift')
        option.rateshift = 10;    
    end    
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrate2(N(1:k),errH1(1:k),option.rateshift,'-*','|| Du-Du_h||',...
              N(1:k),errL2(1:k),option.rateshift,'k-+','|| u-u_h||');
    subplot(1,2,2)
    showrate2(N(1:k),erruIuh(1:k),option.rateshift,'m-+','||Du_I-Du_h||',...
              N(1:k),errMax(1:k),option.rateshift,'r-*','||u_I-u_h||_{\infty}');
end

%% Output
err = struct('N',N(1:k),'H1',errH1(1:k),'L2',errL2(1:k),...
             'uIuhH1',erruIuh(1:k),'uIuhMax',errMax(1:k),'eta',erreta(1:k));
time = struct('N',N(1:k),'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k),...
              'estimate',estimateTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));

%% Display error and CPU time
if option.dispflag
    if ~isfield(option,'dispspace')
       option.dispspace = 5;    % display error for every 5 iterations
    end
    idx = 1:option.dispspace:k;
    
    fprintf('\n');
    disp('Table: Error')
    colname = {'#Dof','||u-u_h||','||Du-Du_h||','||DuI-Du_h||','||uI-u_h||_{max}','eta'};
    disptable(colname,err.N(idx),[],err.L2(idx),'%0.5e',err.H1(idx),'%0.5e',...
          err.uIuhH1(idx),'%0.5e',err.uIuhMax(idx),'%0.5e',err.eta(idx),'%0.5e');
% 
%     disp('Table: CPU time')
%     colname = {'#Dof','Assemble','Solve','Error','Mesh'};
%     disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
%                       time.err,'%0.2e',time.mesh,'%0.2e');     
    fprintf('\n');
end