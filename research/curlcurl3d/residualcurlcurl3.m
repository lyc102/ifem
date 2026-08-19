function [tjump, njump, resK, divJK] = residualcurlcurl3(node,elem,uh,pde)
%% Residual-type error estimator for the 3D Maxwell/H(curl) interface problem
% The discrete solution uh is in the lowest-order ND0 space, using the
% unscaled edge basis
%     phi_k = lambda_i * Dlambda_j - lambda_j * Dlambda_i.
%
% This estimator is based on the decomposition of the error into:
%   - irrotational part: jump of beta u_h and divergence residual,
%   - solenoidal part: element residual of J - curl(mu^{-1}curl u_h) - beta u_h,
%   - tangential jump of mu^{-1} curl u_h.
%
% The corresponding local indicators are of the form
%   eta_0^T = h_T * beta^{-1/2} * ||div J - div(beta u_h)||_{0,T},
%   eta_0^F = h_F^{1/2} * beta_avg^{-1/2} * ||[beta u_h * n]||_{0,F},
%   eta_1^T = h_T * mu^{1/2} * ||J - curl(mu^{-1}curl u_h) - beta u_h||_{0,T},
%   eta_1^F = h_F^{1/2} * mu_avg^{1/2} * ||[n x mu^{-1}curl u_h]||_{0,F}.
%
% For the element residual, note that the Nicaise version differs from the
% Hiptmair formulation by adding face contributions to the element term.
% The implementation below follows the residual form used in the recovery-based
% and interface-adaptive H(curl) analysis.
%
%   References: 
%   Cai, Zhiqiang, and Shuhao Cao. "A recovery-based a posteriori error estimator for H(curl) interface problems." Computer Methods in Applied Mechanics and Engineering 296 (2015): 169-195.
%   Beck, Rudi, Ralf Hiptmair, Ronald HW Hoppe, and Barbara Wohlmuth. "Residual based a posteriori error estimators for eddy current computation." ESAIM: Mathematical Modelling and Numerical Analysis 34, no. 1 (2000): 159-182.
%
%  Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure for edges and faces
NT = size(elem,1);
localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];

% Data structure of the edges
[elem2edge,~,signedge] = dof3edge(elem);

% Data structure of the faces
T = auxstructure3(elem);
face = T.face;
elem2face = T.elem2face;
face2elem = T.face2elem;
neighbor = T.neighbor;
clear T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Dlam,volume] = gradbasis3(node,elem);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NF = size(face,1);
uh2elem = uh(elem2edge).*double(signedge);

[lambda,w] = quadpts3(3); % quadrature order is 3
nQuad = size(lambda,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;

if ~isempty(pde.mu) && isnumeric(pde.mu)
    mu = repmat(pde.mu,NT,1);                
else                            
    mu = pde.mu(center); 
end

if ~isempty(pde.beta) && isnumeric(pde.beta)
    beta = repmat(pde.beta,NT,1);      
else                            
    beta = pde.beta(center);              
end



%%  Data structure for computing the residual
[curluh,~,~] = curlu3(node,elem,uh);

% face2localVertex(k,:,1) is the face k's 
% vertices' local index in face2elem(:,1)
face2localVertex(:,:,1) = localFace(face2elem(:,3),:);
face2localVertex(:,:,2) = localFace(face2elem(:,4),:);

% fn is the global face normal following the index of face matrix
% fn is scaled normal since |fn| = 2*faceArea
fn = cross(node(face(:,2),:) - node(face(:,1),:),...
    node(face(:,3),:) - node(face(:,2),:),2);
% ufn = fn./repmat(sqrt(sum(fn.^2,2)),1,3);
faceArea = 0.5*sqrt(sum(fn.^2,2));
% ufn = fn./repmat(faceArea,1,3)/2;

tic;
%% Compute the element residual and the divergence residual
divJK = zeros(NT,1);
elemresK = zeros(NT,1);

for p = 1:nQuad
    rxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    
    if ~isempty(pde.J) && isnumeric(pde.J)
        Jp = repmat(pde.J,NT,1);
    elseif ~isnumeric(pde.J)
        Jp = pde.J(rxyz);
    else
        Jp = zeros(NT,1);
    end
    
    
    if ~isempty(pde.divJ) && isnumeric(pde.divJ)
        divJKp = repmat(pde.divJ,NT,1);
    elseif ~isnumeric(pde.divJ)
        divJKp = pde.divJ(rxyz);
    else
        divJKp = zeros(NT,1);
    end
    
    
    buhp = zeros(NT,3);
    for k = 1:6
        i = localEdge(k,1); j = localEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        buhp = buhp + ...
            repmat(uh2elem(:,k),1,3).*repmat(beta,1,3).*...
            (lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i));
    end

    elemresK = elemresK + w(p)*sum((Jp - buhp).^2,2);
    divJK = divJK + w(p)*sum(divJKp.^2,2);    
end

divJK = sqrt(divJK./abs(beta)).*volume.^(5/6);

%%% Hiptmair's version
% resK = sqrt(elemresK.*abs(mu)).*volume.^(5/6);

% %%% Nicaise's version
% resK = sqrt(elemresK.*volume.*min(beta.^(-1/2), mu.^(1/2).*volume.^(2/3)));

%%% The element residual that comes from the efficiency proof
resK = sqrt(elemresK.*volume).*...
    sqrt(min(1./abs(beta),abs(mu).*volume.^(2/3)));


%% Compute the norm jump of \beta u_h

% bavg is the avg of the beta across the face following face's index
% see Hiptmair's paper
bavg = 0.5*(beta(face2elem(:,1)) + beta(face2elem(:,2)));

% Compute the normal jump for each face
% Use 2D quadrature points weight
[lambdaFace,wFace] = quadpts(2); % 3 points quad exact for linear vector fields
nQuadFace = size(lambdaFace,1);

njumpFace = zeros(NF,1); 
for p = 1:nQuadFace
    % buh1 is the beta u_h at face2elem(:,1)
    % evaluated at each quadpts on face
    buh1 = zeros(NF,3);
    % buh2 is the beta u_h at face2elem(:,2)
    % evaluated at each quadpts on face
    buh2 = zeros(NF,3);
    for k = 1:6
        i = localEdge(k,1); j = localEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        % lambdaFacei is the reconstructed barycentric from [i1 i2 i3]
        % For example: localEdge = [1 2] and the face we are interested in is
        % locFace = [1 4 3], then lambda(quadp,2) is 0 on this face, data
        % structure is designed to compute face jump for quadrature
        % quantity like the ND0 finite element solution, it circumvents the
        % element-based faceintegral you have to do localFace one by one
        
        % the ending is 1 means the quantity is associated with face2elem(:,1)
        lambdaFacei1 = sum(repmat(lambdaFace(p,:),NF,1)...
            .*(face2localVertex(:,:,1)==i),2);
        lambdaFacej1 = sum(repmat(lambdaFace(p,:),NF,1)...
            .*(face2localVertex(:,:,1)==j),2);
        
        buh1 = buh1 + ...
            repmat(uh2elem(face2elem(:,1),k),1,3).*...
            repmat(beta(face2elem(:,1)),1,3).*...
            (repmat(lambdaFacei1,1,3).*Dlam(face2elem(:,1),:,j)-...
            repmat(lambdaFacej1,1,3).*Dlam(face2elem(:,1),:,i));
        
        lambdaFacei2 = sum(repmat(lambdaFace(p,:),NF,1)...
            .*(face2localVertex(:,:,2)==i),2);
        lambdaFacej2 = sum(repmat(lambdaFace(p,:),NF,1)...
            .*(face2localVertex(:,:,2)==j),2);
        
        buh2 = buh2 + ...
            repmat(uh2elem(face2elem(:,2),k),1,3).*...
            repmat(beta(face2elem(:,2)),1,3).*...
            (repmat(lambdaFacei2,1,3).*Dlam(face2elem(:,2),:,j)-...
            repmat(lambdaFacej2,1,3).*Dlam(face2elem(:,2),:,i));
    end
    % njump_tmp = dot((buh1 - buh2),ufn,2).^2;
    % here use scaled normal fn since |scaled normal| = 2*faceArea
    njump_tmp = 0.25*dot((buh1 - buh2),fn,2).^2;
    njumpFace = njumpFace + wFace(p)*njump_tmp;
end
clear buh1 buh2 njump_tmp lambdaFacei1 lambdaFacei2 ...
    lambdaFacej1 lambdaFacej2
njumpFace = njumpFace./bavg./sqrt(2*faceArea);
njump = 0.5*sum(njumpFace(elem2face),2);
njump = sqrt(njump);


%% Compute the tangential jump of mu^{-1} curl u_h using neighbor


% invmucurluh = mu^{-1} curl u_h
invmucurluh = repmat(1./mu,1,3).*curluh;

% nf vector
normal = -3*repmat(volume,[1,3,4]).*Dlam;

% mu average on each face
invmuavg1 = (1./mu+1./mu(neighbor(:,1),:));
invmuavg2 = (1./mu+1./mu(neighbor(:,2),:));
invmuavg3 = (1./mu+1./mu(neighbor(:,3),:));
invmuavg4 = (1./mu+1./mu(neighbor(:,4),:));

% for boundary edges e, neighbor(t,e) = t. So the difference is zero.     
tjumpFace1 = cross((invmucurluh-invmucurluh(neighbor(:,1),:)),normal(:,:,1),2);
tjumpFace2 = cross((invmucurluh-invmucurluh(neighbor(:,2),:)),normal(:,:,2),2);
tjumpFace3 = cross((invmucurluh-invmucurluh(neighbor(:,3),:)),normal(:,:,3),2);
tjumpFace4 = cross((invmucurluh-invmucurluh(neighbor(:,4),:)),normal(:,:,4),2);

% L2 norm of tangential jump on faces
tjumpFace1 = sum(tjumpFace1.^2,2)./invmuavg1;
tjumpFace2 = sum(tjumpFace2.^2,2)./invmuavg2;
tjumpFace3 = sum(tjumpFace3.^2,2)./invmuavg3;
tjumpFace4 = sum(tjumpFace4.^2,2)./invmuavg4;

tjump = sqrt(tjumpFace1 + tjumpFace2 + tjumpFace3 + tjumpFace4);


%%
recover_time=toc;
fprintf('Time to compute the residual type estimator %6.4g s\n',recover_time);
end