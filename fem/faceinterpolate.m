function uI = faceinterpolate(u,node,elem,elemType)
%% FACEINTERPOLATE interpolate to face elements (RT0 or BDM1).
%
% uI = faceinterpolate(u,node,elem,elemType) interpolates a given function u
% into the lowesr order RT0 or BDM1 finite element spaces. The coefficient
% is given by the line integral int_e u*n ds. The input elemType can be 'RT0'
% or 'BDM1'.
%
% Example: RT0
%  
%      maxIt = 5;
%      node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
%      elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
%      bdFlag = setboundary(node,elem,'Dirichlet');
%      pde = mixBCdata;
%      err = zeros(maxIt,2); 
%      h = zeros(maxIt,1);
%      for k =1:maxIt
%        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%        [u,sigma,eqn] = PoissonRT0(node,elem,bdFlag,pde);
%        sigmaI = faceinterpolate(pde.Du,node,elem,'RT0');
%        err(k,1) = getL2errorRT0(node,elem,pde.Du,sigmaI);
%        err(k,2) = sqrt((sigma-sigmaI)'*eqn.M*(sigma-sigmaI));
%        h(k) = 1./(sqrt(size(node,1))-1);
%      end
%      figure;
%      showrateh2(h,err(:,1),2,'r-+','|| \sigma - \sigma_I ||',...
%                 h,err(:,2),2,'b-+','|| \sigma_h - \sigma_I ||');
%
% Example: BDM1
%
%      maxIt = 5;
%      node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
%      elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
%      bdFlag = setboundary(node,elem,'Dirichlet');
%      pde = mixBCdata;
%      err = zeros(maxIt,2); 
%      h = zeros(maxIt,1);
%      for k =1:maxIt
%        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%        [u,sigma,eqn] = PoissonBDM1(node,elem,bdFlag,pde);
%        sigmaI = faceinterpolate(pde.Du,node,elem,'BDM1');
%        err(k,1) = getL2errorBDM1(node,elem,pde.Du,sigmaI);
%        err(k,2) = sqrt((sigma-sigmaI)'*eqn.M*(sigma-sigmaI));
%        h(k) = 1./(sqrt(size(node,1))-1);
%      end
%      figure;
%      showrateh2(h,err(:,1),2,'r-+','|| \sigma - \sigma_I ||',...
%                 h,err(:,2),2,'b-+','|| \sigma_h - \sigma_I ||');
%
%
% Example: RT1
%
%      maxIt = 5;
%      node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
%      elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
%      bdFlag = setboundary(node,elem,'Dirichlet');
%      pde = mixBCdata;
%      err = zeros(maxIt,2); 
%      h = zeros(maxIt,1);
%      for k = 1:maxIt
%        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%        sigmaI = faceinterpolate(pde.Du,node,elem,'RT1');
%        err(k,1) = getL2errorRT1(node,elem,pde.Du,sigmaI);
%        h(k) = 1./(sqrt(size(node,1))-1);
%      end
%      figure;
%      showrateh2(h,err(:,1),2,'r-+','|| \sigma - \sigma_I ||',...
%                 h,err(:,2),2,'b-+','|| \sigma_h - \sigma_I ||');
%
%
% See also edgeinterpolate, edgeinterpolate1, edgeinterpolate2
%
% Created by Ming Wang at Mar 29, 2011, M-lint modified at May 14, 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('elemType','var'), elemType = 'RT0'; end
%% edge 
if size(elem,2) == 2 % the input elem is edge
    edge = elem;
else
    [elem2edge,edge] = dofedge(elem);
end
NE = size(edge,1);
NT = size(elem,1);
edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
nVec = zeros(NE,2);
nVec(:,1) = edgeVec(:,2); 
nVec(:,2) = -edgeVec(:,1);
clear edgeVec

%% dof for RT0
[lambda,weight] = quadpts1(4);
nQuad = size(lambda,1);
uI = zeros(NE,1);
for i = 1:nQuad
    pxy = lambda(i,1)*node(edge(:,1),:)+lambda(i,2)*node(edge(:,2),:);
    flux = u(pxy);
    uI = uI + weight(i)*dot(flux,nVec,2);   
end

%% dof for BDM1
if strcmp(elemType,'BDM1') || strcmp(elemType,'RT1')
   uI(NE+1:2*NE) = zeros(NE,1);
   for i = 1:nQuad
        pxy = lambda(i,1)*node(edge(:,1),:)+lambda(i,2)*node(edge(:,2),:);
        uI(NE+1:2*NE) = uI(NE+1:2*NE)+ ...
                     weight(i)*3*(lambda(i,1)-lambda(i,2))*dot(u(pxy),nVec,2); 
   end
end

%% dof for RT1
if strcmp(elemType,'RT1')
    mid = (node(edge(:,1),:) + node(edge(:,2),:))/2;
    umid = u(mid);
    uI(2*NE+2*NT) = 0;
    % Face dof coefficients
    uquadpts = umid(elem2edge(:,1),:)+umid(elem2edge(:,2),:)+umid(elem2edge(:,3),:);
    lf = [4*dot(nVec(elem2edge(:,2),:),uquadpts,2) ...
          4*dot(nVec(elem2edge(:,3),:),uquadpts,2)];
    elem2edgeDofValue = [uI(elem2edge(:,1:3)) uI(elem2edge(:,1:3)+NE)];
    localMatrix = [4 8 4 -4  0 4; ...
                   8 4 -4 0 -4 4]';
    lf = lf - elem2edgeDofValue*localMatrix;
    uI((2*NE+1):end) = [2*lf(:,1) - lf(:,2); ...
                        2*lf(:,2) - lf(:,1)]/3;    
end
