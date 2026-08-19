function [errEnergy, errK]= getEnergyerrorNd0boundary(node,elem,pde,uh)
% compute the energy error of the H(curl) problem 
% integration by parts is performed for the curl part of the error
% errK^2=\|\mu^{-1/2} curl(u - u_h)\|_K^2 +  \|\beta^{1/2}(u - u_h)\|_K^2
% = \int_K (f - \beta u_h)\cdot (u - u_h)
% + \int_{\partial K} (mu^{-1} curl(u-u_h)\times n)\cdot (u-u_h) dS


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

%% (I) = \int_K (f - \beta u_h)\cdot (u - u_h)
fbetauhuuh = zeros(NT,1);

[lambda,w] = quadpts3(5); % quadrature order is 5
nQuad = size(lambda,1);


for p = 1:nQuad
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    
    up = pde.exactu(pxyz);
    Jp = pde.J(pxyz);
    uhp = zeros(NT,3);

    for k = 1:6
        i = locEdge(k,1); j = locEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        uhp = uhp + ...
            repmat(uh2elem(:,k),1,3).*...
            (lambda(p,i)*Dlam(:,:,j)-lambda(p,j)*Dlam(:,:,i));
    end
    
    fbetauhuuh = fbetauhuuh + w(p)*...
        dot(up - uhp,Jp - repmat(abs(beta),1,3).*uhp,2);
    
end

fbetauhuuh(isnan(fbetauhuuh)) = 0;
fbetauhuuh = fbetauhuuh.*volume;

%% (II) = + \int_{\partial K} (mu^{-1} curl(u-u_h)\times n)\cdot (u-u_h) dS
invmucurluuhxnuuh = zeros(NT,1);

[lambdaFace,wFace] = quadpts3face(4);
nQuadFace = size(lambdaFace,1);

% normal = |area of the face|*(unit outer normal wrt elem)
normal = -3*repmat(volume,[1,3,4]).*Dlam;

% nQuadFace = # of quad points on each faces
nQuadPerFace = nQuadFace/4;

for p = 1:nQuadFace
    pxyz = lambdaFace(p,1)*node(elem(:,1),:) ...
        + lambdaFace(p,2)*node(elem(:,2),:) ...
        + lambdaFace(p,3)*node(elem(:,3),:) ...
        + lambdaFace(p,4)*node(elem(:,4),:);
    
    if isnumeric(pde.curlu) % a constant vector
        curlup = repmat(pde.curlu,NT,1);
    else
        curlup = pde.curlu(pxyz);
    end
    up = pde.exactu(pxyz);
    
    uhp = zeros(NT,3);
    
    for k = 1:6
        i = locEdge(k,1); j = locEdge(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        uhp = uhp + ...
            repmat(uh2elem(:,k),1,3).*...
            (lambdaFace(p,i)*Dlam(:,:,j)-lambdaFace(p,j)*Dlam(:,:,i));
    end
    invmucurluuhxn = cross(curlup-curluhp,...
        normal(:,:,ceil(p/nQuadPerFace)),2)./repmat(abs(mu),1,3);
    
    invmucurluuhxnuuh = invmucurluuhxnuuh + wFace(p)*...
        dot(invmucurluuhxn,up-uhp,2);
end

invmucurluuhxnuuh(isnan(invmucurluuhxnuuh)) = 0;

%%
% Energy error on each element
errK = sqrt(abs(fbetauhuuh + invmucurluuhxnuuh));

% Total energy norm
% errEnergy =  sqrt(sum(errCurl + errL2));
errEnergy = norm(errK);
% 
% %%%%%%% For solution has strong singularity
% theta = 0.05;
% sortedErr = sort(errCurl + errL2,'descend'); 
% errEnergy = sqrt(sum(sortedErr(floor(theta*NT):end)));
