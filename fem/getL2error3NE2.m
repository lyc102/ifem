function err = getL2error3NE2(node,elem,exactE,Eh,markedElem)
%% GETL2ERROR3NE2 L2 norm of the quadratic (1st type) Nedelec edge element.
% 
%  error = getL23errorNE2(node,elem,exactE,Eh,markedElem);
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Sort elem to ascend ordering
elem = sort(elem,2);

%% Construct Data Structure
[Dlambda,volume] = gradbasis3(node,elem);
% elem2dof
[elem2edge,edge] = dof3edge(elem);
[elem2face,face] = dof3face(elem);
NT = size(elem,1);  NE = size(edge,1); NF = size(face,1);
elem2dof = [elem2edge elem2edge+NE elem2face+2*NE elem2face+2*NE+NF];
% local indices 
locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % phi
               1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % psi
               3 2 4; 3 1 4; 2 1 4; 2 1 3; ...
               4 2 3; 4 1 3; 4 1 2; 3 1 2]; % chi
%% Compute square of the L2 error
err = zeros(NT,1);
[lambda,w] = quadpts3(4); % quadrature order is 4
nQuad = size(lambda,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ... 
        + lambda(p,3)*node(elem(:,3),:) ... 
        + lambda(p,4)*node(elem(:,4),:);
    Ep = exactE(pxy);
    Ehp = zeros(NT,3);
    % compute Ehp at quadrature points
    for k = 1:20
        k1 = locBasesIdx(k,1); 
        k2 = locBasesIdx(k,2); 
        k3 = locBasesIdx(k,3);
        % evaluate basis at quadrature points
        if k<=6
        % phi_k = lambda_{k1}Dlambda_{k2} - lambda_{k2}Dlambda_{k1};
            basis_k = (lambda(p,k1)*Dlambda(:,:,k2) ...
                      -lambda(p,k2)*Dlambda(:,:,k1));
        elseif k<=12
        % phi_k = lambda_{k1}Dlambda_{k2} + lambda_{k2}Dlambda_{k1};
            basis_k = (lambda(p,k1)*Dlambda(:,:,k2) ...
                      +lambda(p,k2)*Dlambda(:,:,k1));
        else
        % chi_k = lambda_{k1}phi_{k2,k3};    
            basis_k = lambda(p,k1)*(lambda(p,k2)*Dlambda(:,:,k3) ...
                                   -lambda(p,k3)*Dlambda(:,:,k2));
        end
        Ehp = Ehp + repmat(Eh(elem2dof(:,k)),1,3).*basis_k;
    end
    err = err + w(p)*volume.*sum((Ep - Ehp).^2,2);
end
% modify the error
err(isnan(err)) = 0; % remove the singular part
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % L2 error on some marked region
end
err = sqrt(sum(err));
%% TODO add example in M-lint