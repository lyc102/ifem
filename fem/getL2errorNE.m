function err = getL2errorNE(node,elem,exactE,Eh,markedElem)
%% GETL2ERRORNE L2 norm of the lowest order Nedelec edge element
% 
%  err = getL2errorNE(node,elem,exactE,Eh) computes the L2 norm between
%  vector field exactE and approximate one from the lowest order Nedelec
%  edge elemnt.
%  
%  err = getL2errorNE(node,elem,exactE,Eh,markedElem) computes the error on
%  marked elements only.
%
% Example
% 
%     [node, elem] = squaremesh([-1 1 -1 1],1);
%     maxIt = 4;
%     L2Err = zeros(maxIt,1);
%     N = zeros(maxIt,1);
%     for k = 1:maxIt
%         [node,elem] = uniformrefine(node,elem);
%         [elem2dof,edge] = dofedge(elem);
%         pde = HodgeLaplacianEdata1;
%         uI = edgeinterpolate(pde.u,node,edge);
%         L2Err(k) = getL2errorNE(node,elem,pde.u,uI);
%         N(k) = length(uI);
%     end
%     r = showrate(N,L2Err,1,'b-+');
%     legend('||u-u_I||',['N^{' num2str(r) '}'],'LOCATION','Best');
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Sort elem to ascend ordering
elem = sort(elem,2);

%% Construct Data Structure
elem2dof = dofedge(elem);
NT = size(elem,1);
[Dlambda,area] = gradbasis(node,elem);
locEdge = [2 3; 1 3; 1 2];

%% Compute square of the L2 error element-wise
[lambda,w] = quadpts(4); % quadrature order is 4
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    if isnumeric(exactE) && all(exactE == 0)   % zero gu
        Ep = 0;
    else
        pxy = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:);
        Ep = exactE(pxy);
    end
    Ehp = zeros(NT,2);
    for k = 1:3 % for each basis
        i = locEdge(k,1); 
        j = locEdge(k,2);
        % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
        Ehp = Ehp + repmat(Eh(elem2dof(:,k)),1,2).*...
                   (lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i));
    end
    err = err + w(p)*sum((Ep - Ehp).^2,2);
end
err = err.*area;
% modify the error
err(isnan(err)) = 0;
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % L2 err on some marked region
end
err = sqrt(sum(err));