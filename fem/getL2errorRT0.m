function [err,elemErr] = getL2errorRT0(node,elem,exactSigma,sigmah,d)
%% GETL2ERRORRT0 L2 norm of RT0 element.
%
%  The input exactSigma can be a function bundle or vector array. 
%  When exactSigma is a vector array, it should be size NT*2, and 
%  exactSigma = \nabla u (or exactSigma= K\nabla u), vector array of size
%  NT by 2.
% 
%  [err,elemErr] = getL2errorRT0(node,elem,exactSigma,sigmah).
%  err gives the L2 norm between exact flux and approximate one from RT0, 
%  elemErr, an NT by 1 array, gives the square of L2 norm between exact 
%  flux and approximate one from RT0 on each element.
%  
%  [err,elemErr] = getL2errorRT0(node,elem,exactSigma,sigmah,d). 
%  If the coefficient 'd' is inlcuded, then the input should be in the way
%  sigmah     = K\nabla u_h  (including K), 
%  exactSigma =  \nabla u    (function bundle, don't including K),    or 
%  exactSigma = K\nabla u    (vector array, including K).
%  err=||K exactSigma - sigma_h||_{K^(-1)} = |||u-u_h||| is the energy norm.
%  elemErr, an NT by 1 array, gives the square of L2 norm between exact 
%  flux and approximate one from RT0 on each element
%
% Example
%
%     maxIt = 5;
%     node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
%     elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
%     bdFlag = setboundary(node,elem,'Dirichlet');
%     pde = mixBCdata;
%     err = zeros(maxIt,1); 
%     h = zeros(maxIt,1);
%     for k = 1:maxIt
%         [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%         [u,sigma] = PoissonRT0(node,elem,bdFlag,pde);
%         err(k) = getL2errorRT0(node,elem,pde.Du,sigma);
%         h(k) = 1./(sqrt(size(node,1))-1);
%     end
%     r1 = showrateh(h,err,2);
%     legend('||\sigma - \sigma_h||',['h^{' num2str(r1) '}'],...
%            'LOCATION','Best');
% 
% See also getHdiverrorRT0, getL2error3RT0.
%
% Created by Ming Wang at Jan 17, 2011, M-lint modified at May 15, 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if exist('d','var') && ~isempty(d)
    if isreal(d)
        K = d;                  % d is an array
    else                        
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:))/3;
    K = d(center);              % d is a function
    end
else
    K = [];
end

%% Construct Data Structure
[elem2dof,~,elem2edgeSign] = dofedge(elem);
NT = size(elem,1);% Ndof = max(elem2dof(:)); %N = size(node,1); 
[Dlambda,area] = gradbasis(node,elem);
locEdge = [2 3; 3 1; 1 2];

%% Compute square of the L2 error element-wise
[lambda,w] = quadpts(3); % quadrature order is 3
nQuad = size(lambda,1);
err = zeros(NT,1);
rotMat = [0 -1; 1 0]; % rotation matrix for computing rotLambda.
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ... 
        + lambda(p,3)*node(elem(:,3),:);
    if ~isnumeric(exactSigma)
        sigmap = exactSigma(pxy);
        if(nargin >=5)&& ~isempty(K) % multiply coeff. if fun bundle flux.
            sigmap = repmat(K,1,2).*sigmap;
        end
    else
        sigmap = exactSigma;         %  K*\nabla u_h.
    end
    sigmahp = zeros(NT,2);
    for k = 1:3 % for each basis
        i = locEdge(k,1); j = locEdge(k,2);
        % phi_k = lambda_iRot_j - lambda_jRot_i;
        sigmahp = sigmahp + repmat(elem2edgeSign(:,k).*sigmah(elem2dof(:,k)),1,2).*...
                   (lambda(p,i)*Dlambda(:,:,j)*rotMat-lambda(p,j)*Dlambda(:,:,i)*rotMat);
    end
    err = err + w(p)*sum((sigmap - sigmahp).^2,2);
end
err = err.*area;               % ||sigma - K\nabla u_h||^2
if(nargin >=5)&& ~isempty(K)   % ||sigma - K\nabla u_h||^2_{K^(-1)}
    err = err./K;
end
elemErr = err;           % ||sigma - K\nabla u_h||^2_{K^(-1)}
% modify the error
err(isnan(err)) = 0;
err = sqrt(sum(err));
