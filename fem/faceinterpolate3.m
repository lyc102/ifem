function uI = faceinterpolate3(u,node,elem,quadOrder)
%% FACEINTERPOLATE3 interpolate to face elements RT0.
%
% uI = faceinterpolate3(u,node,face) interpolates a given function u
% into the lowesr order RT0. The coefficient is given by the face integral 
% int_f u*n ds. 
%
% uI = faceinterpolate3(u,node,elem) when |face| is not given, the function
% will generate the face by |[elem2face,face] = dof3face(elem)|.
%
% uI = faceinterpolate3(u,node,face,6) the last input is the quadrature
% order; see quadpts. 
%
% Example
%   
%   [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
%   maxIt = 3;
%   pde = mixBCdata3;
%   err = zeros(maxIt,1); 
%   h = zeros(maxIt,1);
%   for i = 1:maxIt
%      [node,elem] = uniformrefine3(node,elem);
%      [elem2face,face] = dof3face(elem);
%      uI = faceinterpolate3(pde.Du,node,face);
%      err(i) = getL2error3RT0(node,elem,pde.Du,uI);
%      h(i) = 2^(-i);
%   end
%   figure;
%   showrateh(h,err,2,'-+','|| u - u_I ||');
%
% See also edgeinterpolate, edgeinterpolate1, edgeinterpolate2
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% if ~exist('elemType','var'), elemType = 'RT0'; end

%% Construct Data Structure
if size(elem,2) == 3 % the input elem is face
    face = elem;
else
    elem = sortelem3(elem); 
    [elem2face,face] = dof3face(elem); 
end

%% normal vector of each face
NF = size(face,1);
v12 = node(face(:,2),:) - node(face(:,1),:);
v13 = node(face(:,3),:) - node(face(:,1),:);
nVec = mycross(v12,v13); % |nVec| = 2*area;

%% assemble the right hand side
if ~exist('quadOrder','var'), quadOrder = 3; end
[lambda,weight] = quadpts(quadOrder);
nQuad = size(lambda,1);
uI = zeros(NF,1);
for p = 1:nQuad
    pxyz = lambda(p,1)*node(face(:,1),:) ...
		 + lambda(p,2)*node(face(:,2),:) ...
		 + lambda(p,3)*node(face(:,3),:);
    flux = u(pxyz);
    uI = uI + weight(p)*dot(flux,nVec,2)/2;
end
% if strcmp(elemType,'BDM1')
% end
