function err = getL2error3NE(node,elem,exactE,Eh,markedElem)
%% GETL2ERROR3NE L2 norm of the lowest order Nedelec edge element.
% 
%  err = getL2error3NE(node,elem,exactE,Eh,markedElem)
%
% Example
% 
%     [node, elem] = cubemesh([-1 1 -1 1 -1 1],1);
%     maxIt = 4;
%     L2Err = zeros(maxIt,1);
%     N = zeros(maxIt,1);
%     for k = 1:maxIt
%         [node,elem] = uniformbisect3(node,elem);
%         [elem2dof,edge] = dof3edge(elem);
%         pde = Maxwelldata2;
%         uI = edgeinterpolate(pde.exactu,node,edge);
%         L2Err(k) = getL2error3NE(node,elem,pde.exactu,uI);
%         N(k) = length(uI);
%     end
%     r = showrate(N,L2Err,1,'b-+');
%     legend('||u-u_I||',['N^{' num2str(r) '}'],'LOCATION','Best');
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Sort elem to ascend ordering
elem = sort(elem,2);

%% Construct Data Structure
elem2dof = dof3edge(elem);
NT = size(elem,1);% Ndof = max(elem2dof(:)); %N = size(node,1); 
[Dlambda,volume] = gradbasis3(node,elem);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

%% Compute square of the L2 error element-wise
[lambda,w] = quadpts3(4); % quadrature order is 4
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    if isnumeric(exactE) && all(exactE == 0)   % zero gu
        Ep = 0;
    else
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:) ... 
             + lambda(p,4)*node(elem(:,4),:);
        Ep = exactE(pxyz);
    end
    Ehp = zeros(NT,3);
    for k = 1:6 % for each basis
        i = locEdge(k,1); 
        j = locEdge(k,2);
        % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
        Ehp = Ehp + repmat(Eh(elem2dof(:,k)),1,3).*...
                   (lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i));
    end
    err = err + w(p)*sum((Ep - Ehp).^2,2);
end
err = err.*volume;
% modify the error
err(isnan(err)) = 0;
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % L2 err on some marked region
end
err = sqrt(sum(err));