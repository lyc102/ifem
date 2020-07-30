function M = getmassmat3(node,elem2dof,volume,type,K)
%% GETMASSMAT Get mass matrix of the finite element space
%
% M = GETMASSMAT(node,elem2dof,volume,type,K) get mass matrix of the finite element
% space specified by elemType.  
%
% The type can be: 
% - P1: full mass matrix for P1 element
% - lump: lumped mass matrix for P1 element

N = size(node,1);
NT = length(volume);
Ndof = double(max(elem2dof(:)));

%% Coefficients and default type
if ~exist('type','var'), type = 'P1'; end
if ~exist('K','var'), K = []; end

%% Assembling
n = size(elem2dof,2);
switch n
    case 4  % P1 element
    if strcmp(type,'lump')
        %% Assemble the mass matrix by the mass lumping
        M = accumarray([elem2dof(:,1);elem2dof(:,2);elem2dof(:,3);elem2dof(:,4)],...
                       [volume;volume;volume;volume]/4,[N,1]);    
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
                      node(elem2dof(:,3),:) + node(elem2dof(:,4),:))/4;
            volume = K(center).*volume;
        elseif exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == NT 
            volume = K.*volume;
        end 
        M = sparse(N,N);
        for i = 1:4
            for j = i:4
                ii = double(elem2dof(:,i));
                jj = double(elem2dof(:,j));
                if (j==i)
                    M = M + sparse(ii,jj,volume/10,N,N);
                else
                    M = M + sparse([ii;jj],[jj;ii],[volume/20; volume/20],N,N);                               
                end                    
            end
        end        
    end
    case 10 % P2 element
        %% Assemble the full mass matrix using nodal basis
        % indexing follows Poisson3P2
        [lambda, w] = quadpts3(2);
        nQuad = size(lambda,1);
        ii = zeros(55*NT,1); % 55: # of upper triangular entries
        jj = zeros(55*NT,1); 
        sM = zeros(55*NT,1);
        phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
        phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
        phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
        phi(:,4) = lambda(:,4).*(2*lambda(:,4)-1);
        phi(:,5) = 4*lambda(:,1).*lambda(:,2);
        phi(:,6) = 4*lambda(:,1).*lambda(:,3);
        phi(:,7) = 4*lambda(:,1).*lambda(:,4);
        phi(:,8) = 4*lambda(:,2).*lambda(:,3);
        phi(:,9) = 4*lambda(:,2).*lambda(:,4);
        phi(:,10)= 4*lambda(:,3).*lambda(:,4);
        
        index = 0;
        for i = 1:10
            for j = i:10
                Mij = 0;
                for p = 1:nQuad; Mij = Mij + w(p)*phi(p,i).*phi(p,j); end
                Mij = Mij.*volume;
                ii(index+1:index+NT) = double(elem2dof(:,i));
                jj(index+1:index+NT) = double(elem2dof(:,j));
                sM(index+1:index+NT) = Mij;
                index = index + NT;
            end
        end
        diagIdx = (ii == jj);   upperIdx = ~diagIdx;
        M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
        MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
        M = M + MU + MU';
        
end