function [errQuTotal, errQuK, Qu]= getL2error3ND1d(node,elem,uh,pde)
%% computes \|eps^{1/2}(Q_h u - u_h)\|
% u_h discrete Nedelec linear space
% see also Maxwell1

NT = size(elem,1);
elem2Vhdof = reshape(1:12*NT, 12, NT)';
[Dphi, volume] = gradbasis3(node,elem);

[lambda,w] = quadpts3(3);
nQuad = size(lambda,1);

errQuK = zeros(NT,1);

DiDj = zeros(NT,4,4);
for i = 1:4
    for j = i:4        
        DiDj(:,i,j) = dot(Dphi(:,:,i),Dphi(:,:,j),2);
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end
%% compute Q_u 

%% mass matrix
M = sparse(12*NT,12*NT);
ML = sparse(12*NT,12*NT);
edge_to_vert = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for i = 1:6
    for j = 1:6
        % localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        i1 = edge_to_vert(i,1); i2 = edge_to_vert(i,2);
        j1 = edge_to_vert(j,1); j2 = edge_to_vert(j,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           - (1+(i1==j2))*DiDj(:,i2,j1) ...
                           - (1+(i2==j1))*DiDj(:,i1,j2) ...
                           + (1+(i2==j2))*DiDj(:,i1,j1));
        M = M + sparse(elem2Vhdof(:,i),elem2Vhdof(:,j),Mij,12*NT,12*NT);
        
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           + (1+(i1==j2))*DiDj(:,i2,j1) ...
                           + (1+(i2==j1))*DiDj(:,i1,j2) ...
                           + (1+(i2==j2))*DiDj(:,i1,j1));
        M = M + sparse(elem2Vhdof(:,i+6),elem2Vhdof(:,j+6),Mij,12*NT,12*NT);
        
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           - (1+(i1==j2))*DiDj(:,i2,j1) ...
                           + (1+(i2==j1))*DiDj(:,i1,j2) ...
                           - (1+(i2==j2))*DiDj(:,i1,j1));
        ML = ML + sparse(elem2Vhdof(:,i+6),elem2Vhdof(:,j),Mij,12*NT,12*NT);
    end
end

M = M + ML + ML';

%% rhs (u, phi)
ut = zeros(NT,12);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    up = pde.exactu(pxyz);
    for k = 1:6
        i = edge_to_vert(k,1); j = edge_to_vert(k,2);
        % phi_k = lambda_i Dlambda_j - lambda_j Dlambda_i;
        phi_k = lambda(p,i)*Dphi(:,:,j)-lambda(p,j)*Dphi(:,:,i);
        rhs = dot(phi_k,up,2);
        ut(:,k) = ut(:,k) + w(p)*rhs;
        % psi_k = lambda_i Dlambda_j + lambda_j Dlambda_i;
        psi_k = lambda(p,i)*Dphi(:,:,j)+lambda(p,j)*Dphi(:,:,i);
        rhs = dot(psi_k,up,2);
        ut(:,k+6) = ut(:,k+6) + w(p)*rhs;
    end
end
ut = ut.*repmat(volume,1,12);
U = accumarray(elem2Vhdof(:),ut(:),[12*NT 1]);

%%
Qu = M\U;
%% compute \|eps^{1/2}(Q_h u - u_h)\|
Qu2elem = Qu(elem2Vhdof);
uh2elem = uh(elem2Vhdof);
errQuh = Qu2elem - uh2elem;

Eps2elem = zeros(NT,3,3);  %#ok<PREALL>
% Eps2elem is a NTx3x3 array so that Eps2elem(i,:,:) is the tensor at the
% center of i-th element
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;

Epstmp = arrayfun(@(rowidx) pde.Eps(center(rowidx,:)), ...
    1:size(center,1), 'UniformOutput',0);
Eps2elem = cat(3,Epstmp{:});
Eps2elem = permute(Eps2elem,[3,1,2]); % switch the element idx to the 1st dim


for p = 1:nQuad
    QuminusUhElem = zeros(NT,3);
    for j = 1:12
        je = j; if j > 6; je = je-6; end
        j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);
        
        if j <= 6
            phi_j = lambda(p,j1)*Dphi(:,:,j2)-lambda(p,j2)*Dphi(:,:,j1);
        else
            phi_j = lambda(p,j1)*Dphi(:,:,j2)+lambda(p,j2)*Dphi(:,:,j1);
        end
        QuminusUhElem = QuminusUhElem + repmat(errQuh(:,j),1,3).*phi_j;
    end
    
    EpsQuminusUh = sum(bsxfun(@times, Eps2elem, QuminusUhElem), 2);
    errQuK = errQuK + w(p)*dot(squeeze(EpsQuminusUh), QuminusUhElem, 2);
end

errQuK = sqrt(errQuK.*volume);
errQuTotal = norm(errQuK);

end