function [soln,eqn,err,time,solver] = femStokesHdiv(mesh,pde,option,varargin)
%% FEMSTOKESHDIV solve Stokes equation by H(div) finite element methods
%
%   FEMSTOKESHDIV computes approximations to the Stokes equation on a
%   sequence of meshes obtained by uniform refinement of the input mesh.
% 
% See also femPoisson, femrateStokes
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
    option.elemType = 'RT0-P0';
end

%% Parameters
option = femoption(option);
elemType = option.elemType;
maxIt = option.maxIt;
maxN = option.maxN;
L0 = option.L0;
refType = option.refType;

%% Generate an initial mesh 
for k = 1:L0
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
erruL2 = zeros(maxIt,1); erruIuhH1 = zeros(maxIt,1); erruInf = zeros(maxIt,1);
errpL2 = zeros(maxIt,1); errpIphL2 = zeros(maxIt,1); errpInf = zeros(maxIt,1);
errpL2re = zeros(maxIt,1);
errwL2 = zeros(maxIt,1); errwIwh = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);  h = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
    t = cputime;
    if strcmp(refType,'red')
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
    if isfield(option,'mesh') && strcmp(option.mesh, 'circle')
        bdNode = findboundary(elem,bdFlag);
        temp = sqrt(node(bdNode,1).^2 + node(bdNode,2).^2);
        node(bdNode,:) = node(bdNode,:)./[temp temp];
    %     [node,elem] = optmesh(node,elem);
        node = meshsmoothing(node,elem,5);
    end
    meshTime(k) = cputime - t;
    % solve the equation
    if strcmp(elemType,'RT0-P0') || strcmp(elemType,'RT0')
         % RT0-P0 mixed FEM 
            [soln,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option);
    elseif strcmp(elemType, 'BDM1B-P0') || strcmp(elemType,'BDM1') % (BDM1+bubble)-P0 mixed FEM
            [soln,eqn,info] = StokesBDM1B(node,elem,bdFlag,pde,option);
    end
    % compute error
    t = cputime;
    uh = soln.u;
    ph = soln.p;
    wh = soln.w;
    % ================== error for velocity ==================
    if isfield(pde,'exactu') 
        switch elemType
            case 'RT0-P0'
                erruL2(k) = getL2errorRT0(node,elem,pde.exactu,uh);
                % interpolation
                uI = faceinterpolate(pde.exactu,node,eqn.edge,'RT0');
                ufreeDof = eqn.ufreeDof;
                u0 = uh(ufreeDof);
                uI0 = uI(ufreeDof);
                erruIuhH1(k) = sqrt((u0-uI0)'*eqn.A*(u0-uI0));
                erruInf(k) = max(abs(u0-uI0));
            case 'BDM1B-P0'
                erruL2(k) = getL2errorBDM1(node,elem,pde.exactu,uh);
                % interpolation
                uI = uh;
                uI(1:2*size(eqn.edge,1)) = faceinterpolate(pde.exactu,node,eqn.edge,'BDM1');
                ufreeDof = eqn.ufreeDof;
                u0 = uh(ufreeDof);
                uI0 = uI(ufreeDof);
                erruIuhH1(k) = sqrt((u0-uI0)'*eqn.A*(u0-uI0));
                erruInf(k) = max(abs(u0-uI0));
        end
    end
    % ================== error for pressure ==================
    if isfield(pde,'exactp')
        area = simplexvolume(node,elem);
        errpL2(k) = getL2error(node,elem,pde.exactp,ph);
        pI = Lagrangeinterpolate(pde.exactp,node,elem,'P0');
        pr = recoverP02P1(node,elem,ph,'LA');
        % --------------------------
        errpL2(k) = getL2error(node,elem,pde.exactp,ph);
        errpL2re(k) = getL2error(node,elem,pde.exactp,pr);
        errpIphL2(k) = sqrt(dot((pI-ph).^2,area));
        errpInf(k) = max(abs(ph-pI));
    end
    % ================== error for vort ==================
    if isfield(pde,'exactw')
        switch elemType
            case 'RT0-P0'
                errwL2(k) = getL2error(node,elem,pde.exactw,wh);
                wI = Lagrangeinterpolate(pde.exactw,node,elem,'P1');
                errwIwh(k) = sqrt((wh-wI)'*eqn.Mv*(wh-wI));
            case 'BDM1B-P0'
                errwL2(k) = getL2error(node,elem,pde.exactw,wh);
                wI = Lagrangeinterpolate(pde.exactw,node,elem,'P2',eqn.edge);
                errwIwh(k) = sqrt((wh-wI)'*eqn.Mv*(wh-wI));
        end
    end
    errTime(k) = cputime - t;
    % ================== record ==================
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(uh) + length(ph);
    h(k) = 1./(sqrt(size(node,1))-1);
    if ~isfield(option,'contour')
        option.contour = 0;
    end
    if option.contour % plot contour         
        tricontour(node,elem,wh,10);
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
       if isfield(option,'viewanglew')
           viewanglew = option.viewanglew;
       else
           viewanglew = [-22, 88];
       end
       figure(1); 
       showsolutionRT(node,elem,uh,viewangleu);
       title('Velocity')
       figure(2);
       subplot(1,2,1);
       showsolution(node,elem,ph,viewanglep);
       title('Pressure')
       subplot(1,2,2);
       showsolution(node,elem,wh,viewanglew);
%        showsolution(node,elem,soln.psi,viewanglew);
       title('Vorticity')       
    end
    if N(k) > maxN
        break;
    end
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.80,0.40]);
    subplot(1,3,1)
    showrateh3(h(1:k),erruL2(1:k),2,'k-+','|| u-u_h||',...
               h(1:k),erruIuhH1(1:k),2,'r-*','|| u_I-u_h||_1',...
               h(1:k),erruInf(1:k),2,'b-*','|| u_I-u_h||_{\infty}');
    title('Error of velocity')
    subplot(1,3,2)
    showrateh3(h(1:k),errpL2(1:k),2,'k-+', '|| p - p_h||',...
               h(1:k),errpL2re(1:k),2,'g-+','|| p- p_h^r||',...
               h(1:k),errpIphL2(1:k),2,'r-+','|| p_I - p_h||');
    title('Error of pressure')           
    subplot(1,3,3)
    showrateh2(h(1:k),errwL2(1:k),2,'k-+','|| w - w_h||',...
               h(1:k),errwIwh(1:k),2,'r-+','|| w_I - w_h||');
    title('Error of vorticity')                     
end

% Output
err = struct('h',h(1:k),'N',N,'uL2',erruL2(1:k),'uInf',erruInf,'uIuhH1',erruIuhH1(1:k),...
             'pL2',errpL2(1:k),'pIphL2',errpIphL2(1:k),'pInf',errpInf(1:k),'pL2re',errpL2re(1:k),...
             'wL2',errwL2(1:k),'wIwhL2',errwIwh(1:k));
time = struct('N',N,'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and time
disp('Table: Error')
% error of u
colname = {'#Dof','h','||u_I-u_h||_1','||u-u_h||','||u_I-u_h||_{max}'};
disptable(colname,err.N,[],err.h,'%0.2e',err.uIuhH1,'%0.5e',err.uL2,'%0.5e',err.uInf,'%0.5e');
% error of p
colname = {'#Dof','h','||p_I-p_h||','||p-p_h||','||p_I-p_h||_{max}','||p_I - p^r_h||'};
disptable(colname,err.N,[],err.h,'%0.2e',err.pIphL2,'%0.5e',err.pL2,'%0.5e',err.pInf,'%0.5e',err.pL2re,'%0.5e');
% error of w
colname = {'#Dof','h','||w_I-w_h||','||w-w_h||'};
disptable(colname,err.N,[],err.h,'%0.2e',err.wIwhL2,'%0.5e',err.wL2,'%0.5e');

disp('Table: CPU time')
colname = {'#Dof','Assemble','Solve','Error','Mesh'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                  time.err,'%0.2e',time.mesh,'%0.2e');