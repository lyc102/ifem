function err = getL2error3RT0(node,elem,sigma,sigmah,markedElem)
%% GETL2ERROR3RT0  L2 norm of RT0 element in 3D.
% 
%  The input sigma can now just be function boundle. sigmah is the 
%  flux through faces. 
%  err = getL2error3RT0(node,elem,sigma,sigmah,markedElem)
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
%      uI = faceinterpolate3(pde.Du,node,elem);
%      err(i) = getL2error3RT0(node,elem,pde.Du,uI);
%      h(i) = 2^(-i);
%   end
%   figure;
%   showrateh(h,err,2,'-+','|| u - u_I ||');
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Construct Data Structure
elem = sortelem3(elem); 
elem2face = dof3face(elem); 
% [elem2face,dofSign] = dof3RT0(elem);
NT = size(elem,1);% Ndof = max(elem2dof(:)); %N = size(node,1); 
[Dlambda,volume] = gradbasis3(node,elem);
locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];

%% Compute square of the L2 error element-wise
[lambda,w] = quadpts3(3); % quadrature order is 3
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ... 
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    sigmap = sigma(pxyz);
    sigmahp = zeros(NT,3);
    for l = 1:4 % for each basis
        i = locFace(l,1); j = locFace(l,2); k = locFace(l,3);
        % phi_k = lambda_iRot_j - lambda_jRot_i;
        sigmahp = sigmahp + repmat(sigmah(elem2face(:,l)),1,3)*2.*...
                   (lambda(p,i)*mycross(Dlambda(:,:,j),Dlambda(:,:,k)) + ...
                    lambda(p,j)*mycross(Dlambda(:,:,k),Dlambda(:,:,i)) + ...    
                    lambda(p,k)*mycross(Dlambda(:,:,i),Dlambda(:,:,j)));
    end
    err = err + w(p)*sum((sigmap - sigmahp).^2,2);
end
err = err.*volume;
% modify the error
err(isnan(err)) = 0;
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % L2 err on some marked region
end
err = sqrt(sum(err));