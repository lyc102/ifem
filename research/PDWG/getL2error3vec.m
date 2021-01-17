function [errQuTotal, errQuK, Qu]= getL2error3vec(node,elem,uh,pde)
%% computes \|eps^{1/2}(Q_h u - u_h)\|
% u_h piecwise constant vector represented in edge vector basis

NT = size(elem,1);
elem2intdof = reshape(1:3*NT, 3, NT)';
[~,volume] = gradbasis3(node,elem);

[lambda,w] = quadpts3(4);
nQuad = size(lambda,1);

errQuK = zeros(NT,1); %#ok<PREALL>

%%
loc_e2v = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
elem2ve = zeros(NT,3,6);
for e = 1:6 % six edges
   elem2ve(:,:,e) = node(elem(:,loc_e2v(e,2)),:)-node(elem(:,loc_e2v(e,1)),:);
end

elem2ve = elem2ve./repmat(sqrt(sum(elem2ve.^2,2)),1,3); % normalize for basis

uhElem = zeros(NT,3);

for k = 1:3 % three bases in each K       
    uhElem = uhElem + repmat(uh(elem2intdof(:,k)),1,3).*elem2ve(:,:,k);
end

%% compute Q_u
Qu = zeros(NT,3);
for p = 1:nQuad
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    
    Qup = pde.exactu(pxyz);
    Qu = Qu + w(p)*Qup;
end

%% compute \|eps^{1/2}(Q_h u - u_h)\|

Eps2elem = zeros(NT,3,3);  %#ok<PREALL>
% Eps2elem is a NTx3x3 array so that Eps2elem(i,:,:) is the tensor at the
% center of i-th element
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;

Epstmp = arrayfun(@(rowidx) pde.Eps(center(rowidx,:)), ...
    1:size(center,1), 'UniformOutput',0);
Eps2elem = cat(3,Epstmp{:});
Eps2elem = permute(Eps2elem,[3,1,2]); % switch the element idx to the 1st dim

EpsUhminusQu = sum(bsxfun(@times, Eps2elem, uhElem - Qu), 2);
errQuK =  dot(squeeze(EpsUhminusQu), uhElem - Qu, 2);

errQuK = sqrt(errQuK.*volume);
errQuTotal = norm(errQuK);

%% compute the exactu's interpolation in the edge basis

end