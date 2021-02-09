function [etaK, eqn] = estimatequadcurl(node,elem,soln,pde,option)
%% Residual type estimator for quad curl equation
%    soln: solution of the quad curl equation
%    soln.u: ND_0 approx of curl(curl u) = weak curl of phi
%    soln.phi: CR-P0 approx of -\Delta phi + \nabla p = curl w; div phi = 0
%              size = (NF,3) 
%    soln.w: N_0 approx of curl(curl w) = f (original RHS)
%
%    eta_1(K)^2 = h_K^2 \|f\|_{K}^2 
%          + \sum_F h_F \| [n x curl w_h] \|_{0,F}^2
% 
%    eta_2(K)^2 = h_K^2 \|curl w_h\|_{K}^2 
%          + \sum_F h_F \| [n x grad phi_h] \|_{0,F}^2
%
%    eta_3(K)^2 = h_K^2 \|curl phi_h\|_{K}^2 
%          + \sum_F h_F \| [n x (phi_h - curl u_h)] \|_{0,F}^2
% 
%    Requirement: element array is sorted in the ascending order locally (row-wise)
%       
%    Reference: Error analysis of a decoupled finite element method for quad-curl problems
%               https://arxiv.org/abs/2102.03396
%
%    See also Maxwell1, Stokes3CRP0
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%%
if ~exist('option','var'), option = []; end
if ~isfield(option, 'quadorder'); option.quadorder = 2; end
% whether to include u_h in the computation
% u_h is only for verification purpose
if ~isfield(option, 'includeU'); option.includeU = true; end
% whether to use mesh weighted eta_1
if ~isfield(option, 'scaling'); option.scaling = false; end
if isfield(pde,'J')
    f = pde.J;
elseif isfield(pde,'f')
    f = pde.f;
elseif isfield(pde,'quadcurlu')
    f = pde.quadcurlu;
end


[lambda,weight] = quadpts3(option.quadorder);
nQuad = size(lambda,1);

%% Data structures; geometric quantities
%         [elem,bdFlag] = sortelem3(elem,bdFlag);
elem = sortelem3(elem);
T = auxstructure3(elem);
elem2face = T.elem2face;
face2elem = T.face2elem;
face = T.face;
neighbor = T.neighbor;
clear T;

localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];

NT = size(elem,1); NF = size(face,1);
[Dlambda,volume] = gradbasis3(node,elem);
hK = (6*volume).^(1/3);

nFace = mycross(node(face(:,2),:) - node(face(:,1),:),...
    node(face(:,3),:) - node(face(:,2),:),2);
hF = (sum(nFace.^2,2)).^(1/4);
%% 
curlwh = soln.curlw; phih = soln.phi; uh = soln.u;

%% \sum_F h_F \| [n x curl w_h] \|_{0,F}^2
% scaled exterior normal vector = n_F |F|
normal = -3*repmat(volume,[1,3,4]).*Dlambda;

tJumpcurlwh = zeros(NT,1);
for j = 1:4
    tjump = mycross((curlwh-curlwh(neighbor(:,j),:)),normal(:,:,j),2);
    tJumpcurlwh = tJumpcurlwh + sum(tjump.^2,2);
end
tJumpcurlwh = tJumpcurlwh./hK;

%% h_K^2 \|f\|_{K}^2 
elemResf = zeros(NT,1);

for p = 1:nQuad
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    fp = f(pxyz);
    elemResf = elemResf + weight(p)*sum(fp.^2,2);
end
elemResf = elemResf.*volume.*hK.^2;

%%
% eta1 = sqrt(tJumpcurlwh); % for the test of singular problems
if option.scaling
    eta1 = sqrt(elemResf + tJumpcurlwh).*hK;
%     eta1 = sqrt(elemResf.*hK + tJumpcurlwh.*hK);
else
    eta1 = sqrt(elemResf + tJumpcurlwh);
end
%%  \sum_F h_F \| [n x grad phi_h] \|_{0,F}^2
tJumpgradphih = zeros(NT,3);
for k = 1:3
    gradphihxk = gradu3CR(node,elem,elem2face,phih(:,k),Dlambda); % (grad phi)_{x_k \in {x,y,z}}
    
    for j = 1:4 % 4 faces
        tjump = mycross((gradphihxk-gradphihxk(neighbor(:,j),:)),normal(:,:,j),2);
        tJumpgradphih(:,k) = tJumpgradphih(:,k) + sum(tjump.^2,2);
    end
end
tJumpgradphih = sum(tJumpgradphih,2)./hK;

%% h_K^2 \|curl w_h\|_{K}^2 
elemRescurlwh = sum(curlwh.^2,2).*volume.*hK.^2;

%%
eta2 = sqrt(elemRescurlwh + tJumpgradphih);

%% \sum_F h_F \| [n x (phi_h - curl u_h)] \|_{0,F}^2
curluh = curlu3(node,elem,uh);

phih2elem = zeros(NT,3,4); % 4 faces in the last index
for xj = 1:3 % each component x_j
    phihXj = phih(:,xj);
    phih2elem(:,xj,:) = phihXj(elem2face);
    % Crouzeix-Raviart shape function 1-3\lambda = 1/4 at center
end

% face2localVertex(k,:,1) is the face k's 
% vertices' local indices in face2elem(:,1)-th element
face2localVertex(:,:,1) = localFace(face2elem(:,3),:);
face2localVertex(:,:,2) = localFace(face2elem(:,4),:);

% Compute the tangential jump for (phi_h - curl u_h)
% 3 points quad exact for linear vector fields
[lambdaFace,wFace] = quadpts(2); %
lambdaFace = 1-3*lambdaFace; % change to CR-basis
nQuadFace = size(lambdaFace,1);

tjumpPhihCurluhFace = zeros(NF,1); 
for pf = 1:nQuadFace
    
    phihFaceElem1 = zeros(NF,3); % phi_h - curl u_h at face2elem(:,1) at quad pts
    phihFaceElem2 = zeros(NF,3); % phi_h - curl u_h at face2elem(:,2) at quad pts  
    
    for k = 1:4 
        % iterating on 4 local faces
        % the vectorization is done for faces instead of elements
        % if globally the f-th face is locally the k-th face in an element
        % the tangential jump on f-th face is be computed here
       
        % localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
        ii = localFace(k,1); jj = localFace(k,2); kk = localFace(k,3);
        % phi_h's x component = (phih_x at k-th face)*1 
        %                      +(phih_x at ii-th face)*(1-3\lambda_{ii})
        %                      +(phih_x at jj-th face)*(1-3\lambda_{jj})
        %                      +(phih_x at kk-th face)*(1-3\lambda_{kk})
        % the first term has no contribution to the tangential jump
        %
        % lambdaFaceii is the reconstructed CR-basis for locally k-th face 
        % which has local indices [ii jj kk] 
        
        % is the f-th face locally k-th face in its 1st neighboring element
        isFacekinElem1 = (face2elem(:,3) == k);
        
        % the appended 1 means the quantity is associated with face2elem(:,1)
%         lambdaFaceii1 = sum(repmat(lambdaFace(pf,:),NF,1)...
%             .*(face2localVertex(:,:,1)==ii),2);
%         lambdaFacejj1 = sum(repmat(lambdaFace(pf,:),NF,1)...
%             .*(face2localVertex(:,:,1)==jj),2);
%         lambdaFacekk1 = sum(repmat(lambdaFace(pf,:),NF,1)...
%             .*(face2localVertex(:,:,1)==kk),2);
        
        tmp1 = phih2elem(face2elem(:,1),:,ii)*lambdaFace(pf,1)...
            + phih2elem(face2elem(:,1),:,jj)*lambdaFace(pf,2)...
            + phih2elem(face2elem(:,1),:,kk)*lambdaFace(pf,3);
        
        phihFaceElem1(isFacekinElem1,:) = tmp1(isFacekinElem1,:);
        
        isFacekinElem2 = (face2elem(:,4) == k);
        
%         lambdaFaceii2 = sum(repmat(lambdaFace(pf,:),NF,1)...
%             .*(face2localVertex(:,:,2)==ii),2);
%         lambdaFacejj2 = sum(repmat(lambdaFace(pf,:),NF,1)...
%             .*(face2localVertex(:,:,2)==jj),2);
%         lambdaFacekk2 = sum(repmat(lambdaFace(pf,:),NF,1)...
%             .*(face2localVertex(:,:,2)==kk),2);
        
        tmp2 =  phih2elem(face2elem(:,2),:,ii)*lambdaFace(pf,1)...
                      + phih2elem(face2elem(:,2),:,jj)*lambdaFace(pf,2)...
                      + phih2elem(face2elem(:,2),:,kk)*lambdaFace(pf,3);
        
        phihFaceElem2(isFacekinElem2,:) = tmp2(isFacekinElem2,:);
    end
%     
    phihMinuscurluh1 = phihFaceElem1 - curluh(face2elem(:,1),:);
    phihMinuscurluh2 = phihFaceElem2 - curluh(face2elem(:,2),:);
    
    % here use scaled normal nFace since |scaled normal| = 2*faceArea
    tjump_tmp = mycross((phihMinuscurluh1 - phihMinuscurluh2),nFace,2);
    tjumpPhihCurluhFace = tjumpPhihCurluhFace + 0.25*wFace(pf)*sum(tjump_tmp.^2,2)./hF;
end

tjumpPhihCurluh = 0.5*sum(tjumpPhihCurluhFace(elem2face),2);


%% cheated way just verify the estimator being correct
% phihCenter = zeros(NT,3); % at each elem center
% for j = 1:3 % each component
%     phihj = phih(:,j);
%     phihj2elem = phihj(elem2face);
%     phihCenter(:,j) = sum(phihj2elem,2)/4;
%     % Crouzeix-Raviart shape function 1-3\lambda = 1/4 at center
% end
% % diffPhihCurluh = phihCenter - curluh;
% diffPhihCurluh = phihCenter;
% tjumpPhihCurluh = zeros(NT,1);
% for j = 1:4
%     tmp = mycross((diffPhihCurluh-diffPhihCurluh(neighbor(:,j),:)),normal(:,:,j),2);
%     tjumpPhihCurluh = tjumpPhihCurluh + sum(tmp.^2,2);
% end
% tjumpPhihCurluh = tjumpPhihCurluh./hK;

%% h_K^2 \|curl phi_h\|_{K}^2 
curlphih = curlu3CR(elem2face,phih,Dlambda);
elemResCurlphih = sum(curlphih.^2,2).*volume.*hK.^2;

%% 
% eta3 = sqrt(elemResCurlphih + tjumpPhihCurluh);
eta3 = sqrt(tjumpPhihCurluh);
%%
if option.includeU
    etaK = sqrt(eta1.^2 + eta2.^2 + eta3.^2);
else
    etaK = sqrt(eta1.^2 + eta2.^2);
end
eqn = struct('eta1',eta1,'eta2',eta2,'eta3',eta3, ...
             'elemRescurlwh', sqrt(elemRescurlwh), 'hK', hK);


