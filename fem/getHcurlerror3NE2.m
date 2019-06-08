function err = getHcurlerror3NE2(node,elem,curlE,Eh,markedElem)
%% GETHCURLERROR3NE2 Hcurl norm of approximation error for the quadratic (1st type) Nedelect element in 3-D.
%
% err = getHcurlerror3NE1(node,elem,curlE,Eh,markedElem);
%
% NOTE: it is identical to getHcurlerror3NE since the added basis has no
% contribution to the curl part. 
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
%         [elem2dof,dofSign,edge] = dof3edge(elem);
%         pde = Maxwelldata2;
%         uI = edgeinterpolate1(pde.exactu,node,edge);
%         HcurlErr(k) = getHcurlerror3NE1(node,elem,pde.curlu,uI);
%         N(k) = length(uI);
%     end
%     r = showrate(N,HcurlErr,1,'b-+');
%     legend('||u-u_I||_{curl}',['N^{' num2str(r) '}'],'LOCATION','Best');
%
% See also getHcurlerror3NE1, getHcurlerror3NE2, getL2error3NE
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Sort elem to ascend ordering
elem = sort(elem,2);

%% Construct Data Structure
[Dlambda,volume] = gradbasis3(node,elem);
% elem2dof
[elem2edge,edge] = dof3edge(elem);
% elem2face = dof3RT0(elem);
[elem2face,face] = dof3face(elem);
NT = size(elem,1);  NE = size(edge,1); NF = size(face,1);
elem2dof = [elem2edge elem2edge+NE elem2face+2*NE elem2face+2*NE+NF];
% local indices 
locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % phi
               1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % psi
               3 2 4; 3 1 4; 2 1 4; 2 1 3; ...
               4 2 3; 4 1 3; 4 1 2; 3 1 2]; % chi

%% compute H1 error element-wise using quadrature rule with order quadOrder
err = zeros(NT,1);
[lambda,w] = quadpts3(4);
nQuad = size(lambda,1);
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
    curlEhp = zeros(NT,3);
    % compute Ehp at quadrature points
    for k = 1:20
        k1 = locBasesIdx(k,1); 
        k2 = locBasesIdx(k,2); 
        k3 = locBasesIdx(k,3);
        if k<=6
            % curl phi = 2*Dlambda_i cross Dlambda_j;
            curlBasis_k = 2*cross(Dlambda(:,:,k1),Dlambda(:,:,k2),2);
        elseif k<=12
            % curl psi = 0;
            curlBasis_k = 0;
        else % chi = lambda_{i1}phi_{i2,i3}
            % curl chi =  Dlambda_{i1}cross phi_{i2,i3} + lambda_{i1}curl phi_{i2,i3}
            curlBasis_k = cross(Dlambda(:,:,k1),lambda(p,k2)*Dlambda(:,:,k3) ...
                                   -lambda(p,k3)*Dlambda(:,:,k2),2) ...
                        + lambda(p,k1)*2*cross(Dlambda(:,:,k2),Dlambda(:,:,k3),2);                
        end
        curlEhp = curlEhp + repmat(Eh(elem2dof(:,k)),1,3).*curlBasis_k;
    end
    err = err + w(p)*volume.*sum((curlEp - curlEhp).^2,2);
end
% modify the error
err(isnan(err)) = 0; % remove the singular part
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % L2 error on some marked region
end
err = sqrt(sum(err));
%% TODO write more M-lint