function [err,elemErr] = getL2errorBDM1(node,elem,exactSigma,sigmah,d)
%% GETL2ERRORBDM1  L2 norm of BDM1 element
% 
%  The input exactSigma can be a function bundle or vector array. 
%  When exactSigma is a vector array, it should be size NT*2*3, and 
%  exactSigma(:,:,i) = \nabla u|_i (or exactSigma(:,:,i) = K\nabla u|_i)
%  (i=1,2,3) represent flux at three vertices of a triangle.
% 
%  [err,elemErr] = getL2errorBDM1(node,elem,exactSigma,sigmah).
%  err gives the L2 norm between exact flux and approximate one from BDM1, 
%  elemErr, an NT by 1 array, gives the square of L2 norm between exact 
%  flux and approximate one from BDM1 on each element.
% 
%  [err,elemErr] = getL2errorBDM1(node,elem,exactSigma,sigmah,d). 
%  If the coefficient 'd' is inlcuded, then the input should be in the way
%  sigmah     = K\nabla u_h  (including K), 
%  exactSigma =  \nabla u    (function bundle, don't including K),    or 
%  exactSigma(:,:,i) = K\nabla u|_i (vector array, including K).
%     
%  err=||K exactSigma - sigma_h||_{K^(-1)} = |||u-u_h||| is the energy norm.
%  elemErr, an NT by 1 array, gives the square of L2 norm between exact 
%  flux and approximate one from BDM1 on each element
% 
% Example
%
%     maxIt = 5;
%     node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
%     elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
%     bdEdge = setboundary(node,elem,'Dirichlet');
%     pde = mixBCdata;
%     err = zeros(maxIt,1); N = zeros(maxIt,1);
%     for i =1:maxIt
%         [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
%         [u,sigma,NULL] = PoissonBDM1(node,elem,pde,bdEdge);
%         err(i) = getL2errorBDM1(node,elem,pde.Du,sigma);
%         N(i) = size(u,1);
%     end
%     r1 = showrate(N,err,2);
%     legend('||\sigma - \sigma_h||',['N^{' num2str(r1) '}'],...
%            'LOCATION','Best');
%
% See also getHdiverrorBDM1, getL2errorRT0. 
%
% Created by Ming Wang at Jan 17, 2011, M-lint modified at May 15, 2011. 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 
if nargin >=5 && ~isempty(d)
    if isreal(d)
        K = d;                  % d is an array
    else                        
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:))/3;
    K = d(center);              % d is a function
    end
else
    K=[];
end
%% Construct Data Structure
[elem2dof,edge,elem2edgeSign] = dofedge(elem);
NT = size(elem,1); NE = size(edge,1);
[Dlambda,area] = gradbasis(node,elem);
locEdge = [2 3; 3 1; 1 2];

%% Compute square of the L2 error element-wise
[lambda,w] = quadpts(4); % quadrature order is 3
nQuad = size(lambda,1);
err = zeros(NT,1);
rotMat = [0 -1; 1 0]; % rotation matrix for computing rotLambda.
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ... 
        + lambda(p,3)*node(elem(:,3),:);
    if ~isnumeric(exactSigma)
        sigmap = exactSigma(pxy);
        if(nargin >=5)&& ~isempty(K)
            sigmap = repmat(K,1,2).*sigmap;
        end
    else
        sigmap = lambda(p,1).*exactSigma(:,:,1) + ...
                 lambda(p,2).*exactSigma(:,:,2) + ...
                 lambda(p,3).*exactSigma(:,:,3);
    end
    sigmahp = zeros(NT,2);
    for k = 1:3 % for each basis
        i = locEdge(k,1); j = locEdge(k,2);
        % phi_k = lambda_iRot_j - lambda_jRot_i;
        sigmahp = sigmahp + repmat(elem2edgeSign(:,k).*sigmah(elem2dof(:,k)),1,2).*...
                   (lambda(p,i)*Dlambda(:,:,j)*rotMat-lambda(p,j)*Dlambda(:,:,i)*rotMat)...
                          + repmat(sigmah(NE+elem2dof(:,k)),1,2).*...
                   (lambda(p,i)*Dlambda(:,:,j)*rotMat+lambda(p,j)*Dlambda(:,:,i)*rotMat);     
    end
    err = err + w(p)*sum((sigmap - sigmahp).^2,2);
end
err = err.*area;
if(nargin >=5)&& ~isempty(K)
    err = err./K;
end
elemErr = err;
% modify the error
err(isnan(err)) = 0;
err = sqrt(sum(err));