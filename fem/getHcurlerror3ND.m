function [err, errK] = getHcurlerror3ND(node,elem,curlE,Eh,markedElem)
%% GETHCURLERROR3ND Hcurl norm of approximation error for the lowest order Nedelect element in 3-D.
%
% err = GETHCURLERROR3ND(node,elem,curlE,Eh,markedElem);
%
% Example
% 
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     maxIt = 4;
%     HcurlErr = zeros(maxIt,1);
%     N = zeros(maxIt,1);
%     for k = 1:maxIt
%         [node,elem] = uniformbisect3(node,elem);
%         [elem2dof,edge] = dof3edge(elem);
%         pde = Maxwelldata2;
%         uI = edgeinterpolate(pde.exactu,node,edge);
%         HcurlErr(k) = getHcurlerror3ND(node,elem,pde.curlu,uI);
%         N(k) = length(uI);
%     end
%     r = showrate(N,HcurlErr,1,'b-+');
%     legend('||u-u_I||_{curl}',['N^{' num2str(r) '}'],'LOCATION','Best');
%
% See also getHcurlerror3ND1, getHcurlerror3ND2, getL2error3ND
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Sort elem to ascend ordering
elem = sort(elem,2);

%% Construct Data Structure
elem2dof = dof3edge(elem);
NT = size(elem,1);% Ndof = max(elem2dof(:)); %N = size(node,1); 
% compute curl phi
[Dlambda,volume] = gradbasis3(node,elem);
% curl phi = 2*Dlambda_i cross Dlambda_j;
curlPhi(:,:,6) = 2*cross(Dlambda(:,:,3),Dlambda(:,:,4),2);
curlPhi(:,:,1) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,2),2);
curlPhi(:,:,2) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,3),2);
curlPhi(:,:,3) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,4),2);
curlPhi(:,:,4) = 2*cross(Dlambda(:,:,2),Dlambda(:,:,3),2);
curlPhi(:,:,5) = 2*cross(Dlambda(:,:,2),Dlambda(:,:,4),2);
curlEhp = zeros(NT,3);
for k = 1:6
    curlEhp = curlEhp + ...
    repmat(Eh(elem2dof(:,k)),1,3).*curlPhi(:,:,k);
end

%% compute Hcurl error element-wise
[lambda,w] = quadpts3(2);
nQuad = size(lambda,1);
errK = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    if isnumeric(curlE) % a constant vector
        curlEp = repmat(curlE,NT,1);
    else % function handel
        curlEp = curlE(pxy);
    end
%     curlEp = curlE(pxy);
    % compute Ehp at quadrature points
    errK = errK + w(p)*sum((curlEp - curlEhp).^2,2);
end
errK = errK.*volume;
% modify the error
errK(isnan(errK)) = 0; % remove the singular part
if (nargin == 5) && ~isempty(markedElem)
    errK = errK(markedElem); % error on some marked region
end
err = sqrt(sum(errK));
%% TODO write more M-lint