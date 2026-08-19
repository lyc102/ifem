function [errEnergy, errK]= getEnergyerrorNd0(node,elem,pde,uh)
% compute the energy error of the H(curl) problem
% errK^2=\|\mu^{-1/2} curl(u - u_h)\|_K^2 +  \|\beta^{1/2}(u - u_h)\|_K^2

%% Data structures
[elem2edge,~,signedge] = dof3edge(elem);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
[Dlam,volume] = gradbasis3(node,elem);
NT = size(elem,1); 

%% coefficients
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ND0 basis functions
uh2elem = uh(elem2edge).*double(signedge);
[curluhp,~,~] = curlu3(node,elem,uh);

%% Errors
errL2 = zeros(NT,1);
errCurl = zeros(NT,1);

[lambda,w] = quadpts3(5); % quadrature order is 5
nQuad = size(lambda,1);

for p = 1:nQuad
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    
    up = pde.exactu(pxyz);
    uhp = zeros(NT,3);
    
    if isnumeric(pde.curlu) % a constant vector
        curlup = repmat(pde.curlu,NT,1);
    else 
        curlup = pde.curlu(pxyz);
    end
    
    for k = 1:6
        i = locEdge(k,1); j = locEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        uhp = uhp + ...
            repmat(uh2elem(:,k),1,3).*...
            (lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i));
    end
    
    errCurl = errCurl + w(p)*sum((curlup - curluhp).^2,2)./abs(mu);
    
    errL2 = errL2 + w(p)*sum((up - uhp).^2,2).*abs(beta);
    
end

errL2(isnan(errL2)) = 0;
errL2 = errL2.*volume;
errL2Total2 = sum(errL2);

errCurl(isnan(errCurl)) = 0;
errCurl = errCurl.*volume;

% Energy error on each element
errK = sqrt(errCurl + errL2);

% Total energy norm
% errEnergy =  sqrt(sum(errCurl + errL2));
errEnergy = norm(errK);

if ~all(uh==0)
    fprintf('(L2 error computed by quadrature)/(Total error) is %6g in this refinement. \n',...
        sqrt(errL2Total2)/errEnergy);
else
    fprintf('(L2 norm computed by quadrature)/(Total energy norm) is %6g in this refinement. \n',...
        sqrt(errL2Total2)/errEnergy);
end
% 
% %%%%%%% For solution has strong singularity
% theta = 0.05;
% sortedErr = sort(errCurl + errL2,'descend'); 
% errEnergy = sqrt(sum(sortedErr(floor(theta*NT):end)));