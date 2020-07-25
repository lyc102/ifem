function [soln,eqn,err,time,solver] = femStokes(mesh,pde,option,varargin)
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
% default elemType
if ~isfield(option,'elemType')
    option.elemType = 'CRP0';
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
erruIuhH1 = zeros(maxIt,1); erruL2 = zeros(maxIt,1); erruInf = zeros(maxIt,1);
errpL2 = zeros(maxIt,1); errpIphL2 = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1); h = zeros(maxIt,1);

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
    N(k) = length(uh) + length(ph);
    h(k) = 1./(sqrt(size(node,1))-1);    
    Nv = size(node,1); NT = size(elem,1);    
    % interpolation
    if strcmp(elemType,'P2P0') || strcmp(elemType,'P2P1') || ...
       strcmp(elemType,'isoP2P0') || strcmp(elemType,'isoP2P1')
        uI = pde.exactu([node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2]);
    elseif strcmp(elemType,'CRP0') || strcmp(elemType,'CRP1')
        uI = pde.exactu((node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2);
    elseif strcmp(elemType,'P1bP1')
        % bubble part won't be taken into the error computation
        u0 = pde.exactu(node);
        uI = uh;
        uI([(1:Nv)'; NT+Nv+(1:Nv)']) = u0(:);
    end
    % error for velocity
    erruIuhH1(k) = sqrt((uh-uI(:))'*eqn.Lap*(uh-uI(:)));
    erruInf(k) = max(abs(uh-uI(:)));
    Nu = length(uh);
    uhvec = reshape(uh,Nu/2,2);
    if isfield(pde,'exactu')
        if strcmp(elemType,'P1bP1')
            uhvec = uhvec(1:Nv,:);
        end
        erruL2(k) = getL2error(node,elem,pde.exactu,uhvec);
    end
    % error for pressure
    errpL2(k) = getL2error(node,elem,pde.exactp,ph);
    area = simplexvolume(node,elem);    
    switch elemType(end-1:end)
        case 'P0'
            pI = Lagrangeinterpolate(pde.exactp,node,elem,'P0');
            errpIphL2(k) = sqrt(sum((pI-ph).^2.*area));
        case 'P1'
            pI = Lagrangeinterpolate(pde.exactp,node,elem);
            M = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[Nv,1]);
            ep = pI - ph;
            errpIphL2(k) = sqrt(sum(ep.*(M.*ep)));            
    end
    % record time
    errTime(k) = cputime - t;
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
                errTime(k),erruIuhH1(k), errpL2(k));
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    if ~isfield(option,'contour')
        option.contour = 0;
    end
    if option.contour % plot contour refine it later        
        tricontour(node,elem,uh,10);
    end
    if option.plotflag && length(ph) < 5e4 % show mesh and solution for small size
       if isfield(option,'viewangleu')
           viewangleu = option.viewangleu;
       else
           viewangleu = [25, 87];
       end
       if isfield(option,'viewanglep')
           viewanglep = option.viewanglep;
       else
           viewanglep = [-22, 88];
       end
       figure(1); 
       subplot(1,3,1); showsolution(node,elem,uhvec(:,1),viewangleu);        
       title('Velocity u')
       subplot(1,3,2); showsolution(node,elem,uhvec(:,2),viewangleu);
       title('Velocity v')
       subplot(1,3,3);
       showsolution(node,elem,ph,viewanglep);
       title('Pressure')
    end
    if N(k) > maxN
        break;
    end
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.5,0.40]);
    subplot(1,2,1)
    showrateh3(h(1:k),erruIuhH1(1:k),2,'r-*','| u_I-u_h |_1',...
               h(1:k),erruL2(1:k),2,'k-+','|| u-u_h ||',...
               h(1:k),erruInf(1:k),2,'b-*','|| u_I-u_h ||_{\infty}');
    title('Error of velocity')
    subplot(1,2,2)
    showrateh2(h(1:k),errpL2(1:k),2,'k-+', '|| p - p_h ||',...
               h(1:k),errpIphL2(1:k),2,'r-+','|| p_I - p_h ||');
    title('Error of pressure')           
end

%% Output
err = struct('h',h(1:k),'N',N,'uL2',erruL2(1:k),'uInf',erruInf,'uIuhH1',erruIuhH1(1:k),...
             'pL2',errpL2(1:k),'pIphL2',errpIphL2(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and time
disp('Table: Error')
% error of u
colname = {'#Dof','h','|u_I-u_h|_1','||u-u_h||','||u_I-u_h||_{max}'};
disptable(colname,err.N,[],err.h,'%0.2e',err.uIuhH1,'%0.5e',err.uL2,'%0.5e',err.uInf,'%0.5e');
% error of p
colname = {'#Dof','h','||p_I-p_h||','||p-p_h||'};
disptable(colname,err.N,[],err.h,'%0.2e',err.pIphL2,'%0.5e',err.pL2,'%0.5e');

disp('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');
              