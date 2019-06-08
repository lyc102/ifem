function M = getmassmat(node,elem2dof,area,type,K)
%% GETMASSMAT Get mass matrix of the finite element space
%
% M = GETMASSMAT(elem2edge,area,Dlambda,elemType,K) get mass matrix of the finite element
% space specified by elemType.  
%
% The type can be: 
% - P1: full mass matrix for P1 element
% - lump: lumped mass matrix for P1 element
% - HB: Hierarchical basis for P2 element
% - NB: Nodal basis for P2 element
% - NBB: Nodal basis with bubble functions for P2 element
%
% All cases of P2 elements are added by Lin Zhong.

N = size(node,1);
NT = length(area);
Ndof = double(max(elem2dof(:)));

%% Coefficients and default type
if ~exist('type','var'), type = 'P1'; end
if ~exist('K','var'), K = []; end

%% Assembling
n = size(elem2dof,2);
switch n
    case 3 %% P1 element
    if strcmp(type,'lump')
        %% Assemble the mass matrix by the mass lumping
        M = accumarray([elem2dof(:,1);elem2dof(:,2);elem2dof(:,3)],...
                       [area;area;area]/3,[Ndof,1]);
        if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
            M = K(node).*M;
        elseif exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == N 
            M = K.*M;
        end 
        M = spdiags(M,0,N,N);
    else
        %% Assemble the full mass matrix
        if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
            center = (node(elem2dof(:,1),:) + node(elem2dof(:,2),:) + ...
                      node(elem2dof(:,3),:))/3;
            area = K(center).*area;
        elseif exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == NT 
            area = K.*area;
        end 
        M = sparse(N,N);
        for i = 1:3
            for j = 1:3
                   Mij = area*((i==j)+1)/12;
                   M = M + sparse(elem2dof(:,i),elem2dof(:,j),Mij,Ndof,Ndof);             
            end
        end
    end
    case 6  %% P2 element
    %% Mass matrices for NBB bases
    if strcmp(type,'NBB')
        elem2dof(:,7) = uint32(Ndof+(1:NT)'); % add bubble index    
        Ndof = double(max(elem2dof(:)));
        elemM(1:NT,7)   = 9/20;
        elemM(1:NT,1:3) = 1/20;
        elemM(1:NT,4:6) = 2/15;
        elemM = elemM.*repmat(area,1,7);
        diagM = accumarray(elem2dof(:), elemM(:), [Ndof 1]);
        M = spdiags(double(diagM),0,Ndof,Ndof);
        return;
    end
    %% Mass matrices for HB or NB bases
    [lambda, w] = quadpts(4);
    nQuad = size(lambda,1);
    if strcmp(type,'HB')
        phi(:,1) = lambda(:,1);
        phi(:,2) = lambda(:,2);
        phi(:,3) = lambda(:,3);
    elseif strcmp(type,'NB')
        phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
        phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
        phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    end
    phi(:,4) = 4*lambda(:,2).*lambda(:,3);
    phi(:,5) = 4*lambda(:,3).*lambda(:,1);
    phi(:,6) = 4*lambda(:,1).*lambda(:,2);

    Nbases = 6; NMhalf = 21;
    ii = zeros(NMhalf*NT,1); 
    jj = zeros(NMhalf*NT,1); 
    sM = zeros(NMhalf*NT,1);
    index = 0;
    for i = 1:Nbases
        for j = i:Nbases
            Mij = 0;
            for p = 1:nQuad
                Mij = Mij + w(p)*phi(p,i).*phi(p,j);
            end
            Mij = Mij.*area;
            ii(index+1:index+NT) = double(elem2dof(:,i)); 
            jj(index+1:index+NT) = double(elem2dof(:,j));
            sM(index+1:index+NT) = Mij;
            index = index + NT;
        end
    end
    clear Mij
    diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
    MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
    M = M + MU + MU';
end