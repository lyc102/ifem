function err = getHdiverrorBDM1(node,elem,divSigma,sigmah,markedElem)
%% GETHDIVERRORBDM1 Hdiv norm of the approximation error for BDM1.
%
% The input divSigma is a function bundle, and sigma_h is a vector array 
% whose component including two parts, the first part (i.e., sigmah(1:NE)) 
% is the line integral of flux in norm direction, i.e. 
% sigmah_i = \int_{e_i} sigma\cdot n.
% the second part (i.e., sigmah(NE+1:2*NE)) is the dual basis written as 
% sigmah_i = 3\int_{e_i} (\lambda1-\lambda2) sigma\cdot n  (ei = < 1 , 2 >)
%
% err = getHdiverrorBDM1(node,elem,divSigma,sigmah)
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
%         err(i) = getHdiverrorBDM1(node,elem,pde.f,-sigma,[]);
%         N(i) = size(u,1);
%     end
%     r1 = showrate(N,err,2);
%     legend('||\sigma - \sigma_h||_{H(div)}',['N^{' num2str(r1) '}'],...
%            'LOCATION','Best');
%
% See also getHdiverror3RT0, getL2errorRT0.
%
% Created by Ming Wang at Jan 17, 2011, M-lint modified at May 15, 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct Data Structure
[elem2dof,~,elem2edgeSign] = dofedge(elem);
NT = size(elem,1);
% compute div phi
[Dlambda,area] = gradbasis(node,elem); 
% div phi = 2*(Dlambda_i,Rot_j);
rotMat = [0 -1; 1 0]; % rotation matrix for computing rotLambda.
divPhi(:,3) = 2*dot(Dlambda(:,:,1),Dlambda(:,:,2)*rotMat,2);
divPhi(:,1) = 2*dot(Dlambda(:,:,2),Dlambda(:,:,3)*rotMat,2);
divPhi(:,2) = 2*dot(Dlambda(:,:,3),Dlambda(:,:,1)*rotMat,2);
divSigmahp = zeros(NT,1);
for k = 1:3
    divSigmahp = divSigmahp + elem2edgeSign(:,k).*sigmah(elem2dof(:,k)).*divPhi(:,k);
end

%% compute Hdiv error element-wise
[lambda,w] = quadpts(2);
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    divSigmap = divSigma(pxy);
    % compute divSigmahp at quadrature points
    err = err + w(p)*sum((divSigmap - divSigmahp).^2,2);
end
err = err.*area;
% modify the error
err(isnan(err)) = 0; % remove the singular part
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % error on some marked region
end
err = sqrt(sum(err));
