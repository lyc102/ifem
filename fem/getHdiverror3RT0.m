function err = getHdiverror3RT0(node,elem,divSigma,sigmah,markedElem)
%% GETHDIVERROR3RT0 Hdiv norm of the approximation error for RT0 in 3-D.
%
%  The input divsigma can now just be function boundle. sigmah is the 
%  flux through faces. 
%  err = getHdiverror3RT0(node,elem,divSigma,sigmah) gives the L2 error.
%
% Example
%
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];  % nodes
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7]; % elements
%     elem = label3(node,elem);
%     [node,elem] = uniformbisect3(node,elem);
%     bdFace = setboundary3(node,elem,'Dirichlet','y~=-1','Neumann','y==-1');
%     maxIt = 3;
%     pde = mixBCdata3;
%     err = zeros(maxIt,1); N = zeros(maxIt,1);
%     for i =1:maxIt
%         [node,elem,~,bdFace] = uniformbisect3(node,elem,[],bdFace);
%         barycenter = 1/4.*(node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:)+node(elem(:,4),:));
%         uexa = pde.exactu(barycenter);
%         [u,sigma] = Poisson3RT0(node,elem,pde,bdFace);
%         err(i) = getHdiverror3RT0(node,elem,pde.f,-sigma,[]);
%         N(i) = length(u) + length(sigma) ;
%     end
%     r1 = showrate(N,err,2);
%     legend('|\sigma-\sigma_h|_{div}',['N^{' num2str(r1) '}'],...
%           'LOCATION','Best')
% 
% See also getL2error3RT0.
%
% Created by Ming Wang at Jan 17, 2011, modified m-lint May 15, 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct Data Structure
[elem2dof,dofSign] = dof3RT0(elem);
NT = size(elem,1);
% compute div phi
[Dlambda,volume] = gradbasis3(node,elem); 
% div phi = 3*det(Dlambda_i,Dlambda_j,Dlambda_k);
divPhi(:,4) = 6*dot(Dlambda(:,:,1),cross(Dlambda(:,:,3),Dlambda(:,:,2),2),2);
divPhi(:,1) = 6*dot(Dlambda(:,:,2),cross(Dlambda(:,:,3),Dlambda(:,:,4),2),2);
divPhi(:,2) = 6*dot(Dlambda(:,:,1),cross(Dlambda(:,:,4),Dlambda(:,:,3),2),2);
divPhi(:,3) = 6*dot(Dlambda(:,:,1),cross(Dlambda(:,:,2),Dlambda(:,:,4),2),2);
divSigmahp = zeros(NT,1);
for k = 1:4
    divSigmahp = divSigmahp + ...
    double(dofSign(:,k)).*sigmah(elem2dof(:,k)).*divPhi(:,k);
end

%% compute Hdiv error element-wise
[lambda,w] = quadpts3(2);
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    divSigmap = divSigma(pxy);
    % compute divSigmahp at quadrature points
    err = err + w(p)*sum((divSigmap - divSigmahp).^2,2);
end
err = err.*volume;
% modify the error
err(isnan(err)) = 0; % remove the singular part
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % error on some marked region
end
err = sqrt(sum(err));