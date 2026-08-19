function [eta1, eta2, resK, divJK] = recoverycurlcurl3(node,elem,uh,pde)
% Recovery-type error estimator for Maxwell equation
%   Reference: Cai, Zhiqiang, and Shuhao Cao. "A recovery-based a posteriori error estimator for H(curl) interface problems." Computer Methods in Applied Mechanics and Engineering 296 (2015): 169-195.
%
%  Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure for edges, faces, and neighboring elements
NT = size(elem,1);
localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];

% Data structure of the edges
[elem2edge,edge,signedge] = dof3edge(elem);
elem2edge = double(elem2edge);

% Data structure of the faces and neighbors
T = auxstructure3(elem);
face = T.face;
elem2face = T.elem2face;
clear T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Dlam,volume] = gradbasis3(node,elem);
[Dlam,volume] = gradbasis3(node,elem);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NE = size(edge,1); NF = size(face,1);
uh2elem = uh(elem2edge).*double(signedge);

[curluhp,~,~] = curlu3(node,elem,uh);

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

%%  Data structure for the recovery G(\mu^{-1} curl u_h)


%% M = Mass matrix for ND0 element
M = sparse(NE,NE);
tic;
for i = 1:6
    for j = 1:6
        signij = signedge(:,i).*signedge(:,j);
        ind = (signij == -1);
        
        % localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        i1 = localEdge(i,1); i2 = localEdge(i,2);
        j1 = localEdge(j,1); j2 = localEdge(j,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        Mij = 1/20*volume.*mu.*...
            ((1+(i1==j1))*dot(Dlam(:,:,i2),Dlam(:,:,j2),2)...
            - (1+(i1==j2))*dot(Dlam(:,:,i2),Dlam(:,:,j1),2)...
            - (1+(i2==j1))*dot(Dlam(:,:,i1),Dlam(:,:,j2),2)...
            + (1+(i2==j2))*dot(Dlam(:,:,i1),Dlam(:,:,j1),2));
        Mij(ind) = -Mij(ind);
        M = M + sparse(elem2edge(:,i),elem2edge(:,j),Mij,NE,NE);
    end
end
clear signij ind

gt = zeros(NT,6);
for p = 1:nQuad
    %   localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    for k = 1:6
        i = localEdge(k,1); j = localEdge(k,2);
        phi_k = lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i);
        rrhs = double(signedge(:,k)).*dot(phi_k,curluhp,2);
        gt(:,k) = gt(:,k) + w(p)*rrhs;
    end
end
gt = gt.*repmat(volume,1,6);
G = accumarray(elem2edge(:),gt(:),[NE 1]);
clear gt rrhs


%% L^2 Recovery of sigma = G(\mu^(-1) curl u_h)
sigma = M\G;
% sigma = gmres(M,G,8,1e-7,20);
sigma2elem = sigma(elem2edge).*double(signedge);

%% Compute the Solenoidal error estimator
eta1 = zeros(NT,1);

for p = 1:nQuad
    sigmahp = zeros(NT,3);
    for k = 1:6
        i = localEdge(k,1); j = localEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i.
        
        sigmahp = sigmahp + ...
            repmat(sigma2elem(:,k),1,3).*...
            (lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i));
    end
    
    eta1 = eta1 + w(p)*mu.*sum((sigmahp - curluhp./repmat(mu,1,3)).^2,2);
end
clear sigmahp curluhp

eta1 = sqrt(eta1.*volume);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Face data structure for the recovery G(\beta u_h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e2fArea(:,1) = 3*sqrt(sum(Dlam(:,:,1).^2,2)).*volume;
e2fArea(:,2) = 3*sqrt(sum(Dlam(:,:,2).^2,2)).*volume;
e2fArea(:,3) = 3*sqrt(sum(Dlam(:,:,3).^2,2)).*volume;
e2fArea(:,4) = 3*sqrt(sum(Dlam(:,:,4).^2,2)).*volume;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% another type of signface %%%%%%%%%%%%%%%%%
%%% if the local outer normal coincide with the global
%%% normal of each face, then the sign is 1
facenormal = cross(node(face(:,2),:) - node(face(:,1),:),...
    node(face(:,3),:) - node(face(:,2),:),2);
%n2elem is the global face normal to each element
n2elem(:,:,1) = facenormal(elem2face(:,1),:);
n2elem(:,:,2) = facenormal(elem2face(:,2),:);
n2elem(:,:,3) = facenormal(elem2face(:,3),:);
n2elem(:,:,4) = facenormal(elem2face(:,4),:);
% on2elem is the outward normal on each element
on2elem = -Dlam;
signface = -ones(NT,4,'int8');
signface(:,1) = signface(:,1)+2*int8(dot(n2elem(:,:,1),on2elem(:,:,1),2)>0);
signface(:,2) = signface(:,2)+2*int8(dot(n2elem(:,:,2),on2elem(:,:,2),2)>0);
signface(:,3) = signface(:,3)+2*int8(dot(n2elem(:,:,3),on2elem(:,:,3),2)>0);
signface(:,4) = signface(:,4)+2*int8(dot(n2elem(:,:,4),on2elem(:,:,4),2)>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% R = Mass matrix for RT0 element
R = sparse(NF,NF);
for i = 1:4
    for j = i:4
        % local to global index map and its sign
        ii = double(elem2face(:,i));
        jj = double(elem2face(:,j));
        signij = signface(:,i).*signface(:,j);
        ind = (signij == -1);
        i1 = localFace(i,1); i2 = localFace(i,2); i3 = localFace(i,3);
        % [i1,i2,i3] is the face opposite to vertex i.
        % psi_i = |F_i|/6|T|(lambda_i1 e_i1 + lambda_i2 e_i2 + lambda_i3 e_i3);
        
        j1 = localFace(j,1); j2 = localFace(j,2); j3 = localFace(j,3);
        
        Rij = 1/20*volume./beta.*e2fArea(:,i).*e2fArea(:,j).*( ...
            (1+(i1==j1))*dot(cross(Dlam(:,:,i2),Dlam(:,:,i3),2), ...
            cross(Dlam(:,:,j2),Dlam(:,:,j3),2),2)...
            +(1+(i1==j2))*dot(cross(Dlam(:,:,i2),Dlam(:,:,i3),2), ...
            cross(Dlam(:,:,j3),Dlam(:,:,j1),2),2)...
            +(1+(i1==j3))*dot(cross(Dlam(:,:,i2),Dlam(:,:,i3),2), ...
            cross(Dlam(:,:,j1),Dlam(:,:,j2),2),2)...
            +(1+(i2==j1))*dot(cross(Dlam(:,:,i3),Dlam(:,:,i1),2), ...
            cross(Dlam(:,:,j2),Dlam(:,:,j3),2),2)...
            +(1+(i2==j2))*dot(cross(Dlam(:,:,i3),Dlam(:,:,i1),2), ...
            cross(Dlam(:,:,j3),Dlam(:,:,j1),2),2)...
            +(1+(i2==j3))*dot(cross(Dlam(:,:,i3),Dlam(:,:,i1),2), ...
            cross(Dlam(:,:,j1),Dlam(:,:,j2),2),2)...
            +(1+(i3==j1))*dot(cross(Dlam(:,:,i1),Dlam(:,:,i2),2), ...
            cross(Dlam(:,:,j2),Dlam(:,:,j3),2),2)...
            +(1+(i3==j2))*dot(cross(Dlam(:,:,i1),Dlam(:,:,i2),2), ...
            cross(Dlam(:,:,j3),Dlam(:,:,j1),2),2)...
            +(1+(i3==j3))*dot(cross(Dlam(:,:,i1),Dlam(:,:,i2),2), ...
            cross(Dlam(:,:,j1),Dlam(:,:,j2),2),2));
        Rij(ind) = -Rij(ind);
        
        if (j==i)
            R = R + sparse(ii,jj,Rij,NF,NF);
        else
            R = R + sparse([ii;jj],[jj;ii],[Rij; Rij],NF,NF);
        end
    end
end
clear signij ind

%% Compute (u_h, psi_k)
kt = zeros(NT,4);
for p = 1:nQuad
    uhp = zeros(NT,3);
    
    for m = 1:6
        j1 = localEdge(m,1); j2 = localEdge(m,2);
        uhp = uhp + ...
            repmat(uh2elem(:,m),1,3).*...
            (lambda(p,j1)*Dlam(:,:,j2)-lambda(p,j2)*Dlam(:,:,j1));
    end
    
    %   localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
    for k = 1:4
        i1 = localFace(k,1); i2 = localFace(k,2); i3 = localFace(k,3);
        % psi_k = |F_k|/6|T|(lambda_i1 e_i1elemresK = zeros(NT,1); + lambda_i2 e_i2 + lambda_i3 e_i3);
        psi_k = lambda(p,i1)*cross(Dlam(:,:,i2),Dlam(:,:,i3),2)...
            + lambda(p,i2)*cross(Dlam(:,:,i3),Dlam(:,:,i1),2)...
            + lambda(p,i3)*cross(Dlam(:,:,i1),Dlam(:,:,i2),2);
        
        rrhs = double(signface(:,k)).*e2fArea(:,k).*dot(psi_k,uhp,2);
        kt(:,k) = kt(:,k) + w(p)*rrhs;
    end
end
kt = kt.*repmat(volume,1,4);
K = accumarray(elem2face(:),kt(:),[NF 1]);
clear kt rrhs psi_k


%%  L^2 Recovery of tau = G(\beta u_h)
tau = R\K;
% tau = gmres(R,K,8,1e-7,20);
tau2elem = tau(elem2face).*double(signface);

% [lambda,w] = quadpts3(3); % quadrature order is 3
% nQuad = size(lambda,1);

%% Compute the Irrotational error estimator
eta2 = zeros(NT,1);

for p = 1:nQuad
    uhp = zeros(NT,3);
    
    %   localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
    for m = 1:6
        j1 = localEdge(m,1); j2 = localEdge(m,2);
        uhp = uhp + ...
            repmat(uh2elem(:,m),1,3).*...
            (lambda(p,j1)*Dlam(:,:,j2)-lambda(p,j2)*Dlam(:,:,j1));
    end
    
    tauhp = zeros(NT,3);
    for k = 1:4
        i1 = localFace(k,1); i2 = localFace(k,2); i3 = localFace(k,3);
        % psi_k = |F_k|/6|T|(lambda_k1 e_k1 + lambda_k2 e_k2 + lambda_k3 e_k3);
        
        tauhp = tauhp + ...
            repmat(tau2elem(:,k),1,3).*repmat(e2fArea(:,k),1,3).*...
            (lambda(p,i1)*cross(Dlam(:,:,i2),Dlam(:,:,i3),2)...
            + lambda(p,i2)*cross(Dlam(:,:,i3),Dlam(:,:,i1),2)...
            + lambda(p,i3)*cross(Dlam(:,:,i1),Dlam(:,:,i2),2));
    end
    
    eta2 = eta2 + w(p)*sum((tauhp - repmat(beta,1,3).*uhp).^2,2)./abs(beta);
end
clear uhp tauhp

eta2 = eta2.*volume;

eta2 = sqrt(eta2);

%%  Compute the Recover type element residual(see the paper)
[curlsigma,~,~] = curlu3(node,elem,sigma);

elemresK = zeros(NT,1);
divJK = zeros(NT,1);
% avgRK is the average 1/|K| \int_K (f - \beta u_h - curl \sigma)
avgRK = zeros(NT,3);



% [lambda,w] = quadpts3(5); % quadrature order is 5
% nQuad = size(lambda,1);

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
            repmat(uh2elem(:,k),1,3).*...
            repmat(beta,1,3).*...
            (lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i));
    end
    
    residualK = Jp - buhp - curlsigma;
    avgRK = avgRK + w(p)*residualK;
%     elemresK = elemresK + w(p)*sum((Jp - tauhp - curlsigma).^2,2);
    elemresK = elemresK + w(p)*sum(residualK.^2,2);
    divJK = divJK + w(p)*sum(divJKp.^2,2);
   
    
end
clear Jp divJKp tauhp buhp tau2elem sigma2elem rxyz

%% Subtract the average vector field
elemresavgK = zeros(NT,1);
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
    
    buhp = zeros(NT,3);
    for k = 1:6
        i = localEdge(k,1); j = localEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        buhp = buhp + ...
            repmat(uh2elem(:,k),1,3).*...
            repmat(beta,1,3).*...
            (lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i));
    end
    
    elemresavgK = elemresavgK + w(p)*sum((Jp - buhp - curlsigma - avgRK).^2,2);
end
clear Jp buhp curlsigma avgRK rxyz uh2elem

%% Compute the element residual estimator

%%% If use the averaged version of element residual
% elemresK = elemresavgK;

%%% From the LS point of view
% resK = sqrt(elemresK./abs(beta).*volume);

%%% The element residual that comes from Helmholtz decomposition
% resK = sqrt(elemresK.*abs(mu).*volume.^(1+2/3));

%%% The element residual that comes from the efficiency proof
resK = sqrt(elemresK.*volume).*...
    sqrt(min(1./abs(beta),abs(mu).*volume.^(2/3)));
% resK = sqrt(elemresK.*volume).*...
%     (1./abs(beta) + abs(mu).*volume.^(2/3));

%%% The element residual that comes from flux recovery


divJK = divJK.*volume;
divJK = sqrt(divJK./abs(beta).*volume.^(2/3));

%%
recover_time=toc;
fprintf('Time to recover the error %6.4g s\n',recover_time);
end