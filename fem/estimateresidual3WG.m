function eta = estimateresidual3WG(node,elem,Du,pde)

% Du: graident of CR part

NT = size(elem,1);
%% Geometric quantity
face = uint32([elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])]);
v12 = node(face(:,2),:)-node(face(:,1),:);
v13 = node(face(:,3),:)-node(face(:,1),:);
normal = mycross(v12,v13,2);
normal = reshape(normal,NT,3,4);
v12 = v12(3*NT+1:4*NT,:); 
v13 = v13(3*NT+1:4*NT,:);
v14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = abs(dot(mycross(v12,v13,2),v14,2)/6);
h = volume.^(1/3);
%% Jump of tangential component of gradient
% data structure
T = auxstructure3(elem);
neighbor = T.neighbor; 
clear T
% face jump of Dub
faceJumpb = zeros(NT,1);
for i = 1:4
    faceJumpb = faceJumpb + sum(mycross((Du-Du(neighbor(:,i),:)),normal(:,:,i),2).^2,2); ...
end
faceJumpb = faceJumpb./h; % the normal is h^2
faceJumpo = 0;    % skip the interior part since it is of high order  
%% Elementwise residual
elemResidual = zeros(NT,1);
if isfield(pde,'f') && isnumeric(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;    
    fc = pde.f(center);
    [lambda,weight] = quadpts3(3);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ...
             + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
        fp = pde.f(pxyz);
        elemResidual = elemResidual + weight(p)*(fp-fc).^2;
    end
    elemResidual = elemResidual.*h.^5; % h^2\int _T ||f-fc||^2
end
%% Residual type error estimator
eta = (elemResidual + 0.5*faceJumpb + 0.5*faceJumpo).^(1/2);